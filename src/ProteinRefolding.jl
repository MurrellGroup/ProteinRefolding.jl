module ProteinRefolding

include("Process.jl")
using .ProcessPDB

using ColabMPNN, PyBoltz, PyBoltz.Schema, Combinatorics
import BioStructures
using BioStructures: MolecularStructure
using ProteinChains: ProteinStructure, ProteinChain, readpdb
using Combinatorics: Permutations, nthperm
using TMscore: run_tmscore, tmscore

export refold

#=
TODO
- showprogress
- Fix processed_proteins structure_names
- Allow usage of custom sequence
=#

torch = PyBoltz.PythonCall.pyimport("torch")
torch.torch.set_float32_matmul_precision("high")

function mpnn_generate_seqs(
    pdb_filename::AbstractString, chainids=BioStructures.chainids(readfile(pdb_filename));
    msa_size=1,
    query_temperature=0.1,
    extra_temperature=0.1,
    mpnn_model_args=(),
    mpnn_model_kwargs=(;),
    mpnn_model=ColabMPNN.mk_mpnn_model(mpnn_model_args...; mpnn_model_kwargs...)
)
    ColabMPNN.prep_inputs(mpnn_model, pdb_filename; chain=join(collect(chainids), ','))
    query_seq = ColabMPNN.sample(mpnn_model; temperature=query_temperature).seq |> only
    extra_seqs = ColabMPNN.sample_parallel(mpnn_model, msa_size-1; temperature=extra_temperature).seq

    return query_seq, extra_seqs
end

function make_molecularinput(query_seq, extra_seqs, input)
    msa = [query_seq; extra_seqs]
    query_chains = split(query_seq, "/")
    msa_chains = split.(msa, "/")
    chainids = BioStructures.chainids(readfile(input.pdb_filename))

    molecularinput = Schema.MolecularInput(
        name=first(splitext(basename(input.pdb_filename))),
        sequences = [Schema.protein(id=chainids[ch_i], sequence=query_chains[ch_i],  msa=map(multichainseq -> multichainseq[ch_i], msa_chains)) for ch_i in eachindex(chainids)] 
    )
    return molecularinput
end

function chainwise_run_tmscore(A::ProteinStructure, B::ProteinStructure)
    tm_results = map(A.chains, B.chains) do A_ch, B_ch
        A_ch.id == B_ch.id || @warn "Chain IDs do not match: $A_ch and $B_ch."
        A_struc = MolecularStructure(ProteinStructure("$(A.name)_$(A_ch.id)", [A_ch]))
        B_struc = MolecularStructure(ProteinStructure("$(B.name)_$(B_ch.id)", [B_ch]))
        run_tmscore(A_struc, B_struc)
    end
    return tm_results
end

concat_chains(chains::Vector{<:ProteinChain}, id="A") = ProteinChain(id, reduce(vcat, c.atoms for c in chains), join("A"^length(c) for c in chains))

function complex_run_tmscore!(A::ProteinStructure, B::ProteinStructure; max_seconds=60.0)
    isoptimal = true

    start_time = time()
    B_struc = MolecularStructure(ProteinStructure(B.name, [concat_chains(B.chains)]))
    tmresults = []
    for A_chain_perm in permutations(A.chains)
        length.(A_chain_perm) == length.(B.chains) || continue
        A_struc = MolecularStructure(ProteinStructure(A.name, [concat_chains(A_chain_perm)]))
        push!(tmresults, run_tmscore(A_struc, B_struc))

        if time() - start_time > max_seconds
            @warn "Took more time than max_seconds(=$(max_seconds)). Not running more permutations. Answer will be suboptimal. A: $(A.name). B: $(B.name). nperms: $(length(permutations(A.chains))) "
            isoptimal = false
            break
        end
    end
    n = argmax(i -> tmresults[i].tmscore, eachindex(tmresults))
    nthperm!(A, n)
    return (; tmresult=tmresults[n], isoptimal)
end
complex_run_tmscore(A, B; kws...) = complex_run_tmscore!(deepcopy(A), B; kws...)

function preprocess_and_generate(proteinfiles::Vector{String}, temp_dir; warn_pdbformat=true, seqgen_kws...)
    preprocessed_pdbs = preprocess.(proteinfiles, joinpath(temp_dir, "preprocessed_pdbs"))
    seqs = mpnn_generate_seqs.(getproperty.(preprocessed_pdbs, :pdb_filename); seqgen_kws...)
    molecularinput = make_molecularinput.(first.(seqs), last.(seqs), preprocessed_pdbs)
    fillgaps!.(molecularinput, preprocessed_pdbs; warn_pdbformat)
    return preprocessed_pdbs, molecularinput
end

function predict_and_postprocess(preprocessed_pdbs::Vector{<:PreprocessedPDB}, molecularinput::Vector{MolecularInput}; run_tmscore_func = complex_run_tmscore, devices=1)
    predicted_proteins = predict(molecularinput, ProteinStructure; devices)
    names_in_predicted = Set(getproperty.(predicted_proteins, :name))
    preprocessed_pdbs = filter(pdb -> readfile(pdb.pdb_filename).name ∈ names_in_predicted, preprocessed_pdbs)

    resetgaps!.(predicted_proteins, preprocessed_pdbs; warn_pdbformat=false)
    resetchainids!.(predicted_proteins, preprocessed_pdbs)
    preprocessed_proteins = resetprocessedpdb.(preprocessed_pdbs)

    tmresults = run_tmscore_func.(preprocessed_proteins, predicted_proteins)
    return (; tmresults, preprocessed_proteins, predicted_proteins)
end

# TODO Check that this works in all special cases, like missing gaps etc
function insert_origseq!(molinput::MolecularInput, preprocessed_pdb::PreprocessedPDB)
    prot = readpdb(preprocessed_pdb.pdb_filename)
    for i in eachindex(prot.chains)
        # TODO this assertion shouldn't be needed but I'm keeping it for now
        @assert molinput["sequences"][i]["protein"]["id"] == prot.chains[i].id "$(molinput["sequences"][i]["protein"]["id"]) $(prot.chains[i].id)"
        molinput["sequences"][i]["protein"]["sequence"] = prot.chains[i].sequence
    end
end

function refold(proteinfiles::Vector{String}, temp_dir; reuse_queryseq=false, warn_pdbformat=true, run_tmscore_func = complex_run_tmscore, devices=1, seqgen_kws...)
    preprocessed_pdbs, molecularinput = preprocess_and_generate(proteinfiles, temp_dir; warn_pdbformat, seqgen_kws...)
    reuse_queryseq && insert_origseq!.(molecularinput, preprocessed_pdbs)
    return predict_and_postprocess(preprocessed_pdbs, molecularinput; run_tmscore_func, devices)
end
refold(proteinfiles::Vector{String}; kws...) = mktempdir(temp_dir -> refold(proteinfiles, temp_dir; kws...))
refold(proteinfile::String, args...; kws...) = refold([proteinfile], args...; kws...)

end
