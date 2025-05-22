module ProcessPDB

using BioStructures, BioSequences
using BioStructures: threeletter_to_aa, AA_G, AA_Gap, AminoAcid
using PyBoltz.Schema: MolecularInput
using ProteinChains: ProteinStructure, ProteinChain, readpdb, writecif
import ProteinChains

export 
    preprocess, 
    fillgaps!, 
    resetgaps!, 
    resetchainids!, 
    resetprocessedpdb,
    PreprocessedPDB,
    readfile

struct PreprocessedPDB{T}
    pdb_filename::String
    origchainids::Dict{String, String}
    chainseqs::Dict{String, String}
    chainkeepindices::Dict{String, Vector{Int}}
    origformat::Type{T}
end

fileformat(proteinfile) = ifelse(endswith(proteinfile, ".pdb"), PDBFormat, MMCIFFormat)
readfile(proteinfile) = read(proteinfile, fileformat(proteinfile); structure_name = first(splitext(basename(proteinfile))), remove_disorder=true, read_het_atoms=false)

function convertchainids!(struc)
    origchainids = Dict{String, String}()
    next_chainid = 'A'
    for ch in collect(defaultmodel(struc))
        if isempty(collectatoms(ch, backboneselector))
            continue
        end
        while string(next_chainid) ∈ chainids(struc) 
            next_chainid += 1
        end

        origchainids[string(next_chainid)] = chainid(ch)
        chainid!(ch, string(next_chainid))
    end
    return origchainids
end

function chainresnumbers(struc::MolecularStructure, origchainids)
    chainresnums = Dict{String, Set{Int}}()
    for ch in defaultmodel(struc)
        chainresnums[origchainids[chainid(ch)]] = Set(resnumber.(collectresidues(ch)))
    end
    return chainresnums
end


function find_gapinfo(proteinfile, chainresnums::Dict{String, Set{Int}}, ::Type{MMCIFFormat})
    mmcif = MMCIFDict(proteinfile)

    # Find sequences
    chainseqs = Dict{String, String}()
    for (_chainids, _seq) in zip(mmcif["_entity_poly.pdbx_strand_id"], mmcif["_entity_poly.pdbx_seq_one_letter_code_can"])
        for _chainid in split(_chainids, ",")
            _seq = join(split(_seq))
            chainseqs[_chainid] = _seq
        end
    end

    # Find index mapping
    chainkeepindices = Dict{String, Vector{Int}}()
    _chainids = mmcif["_atom_site.auth_asym_id"]
    _seqindices = tryparse.(Int, mmcif["_atom_site.label_seq_id"])
    _resindices = tryparse.(Int, mmcif["_atom_site.auth_seq_id"])

    for (_chainid, seqind, resind) in zip(_chainids, _seqindices, _resindices)
        (_chainid in keys(chainresnums) && resind in chainresnums[_chainid]) || continue
        keepindices = get!(chainkeepindices, _chainid, Int[])
        push!(keepindices, seqind)
    end
    chainkeepindices = Dict{String, Vector{Int}}(k => sort(unique(v)) for (k, v) in pairs(chainkeepindices))

    return chainseqs, chainkeepindices
end
find_gapinfo(proteinfile, chainresnums) = find_gapinfo(proteinfile, chainresnums, fileformat(proteinfile))

function fillgaps!(molecularinput::MolecularInput, pdb::PreprocessedPDB{MMCIFFormat}; kws...)
    for ch_input in molecularinput["sequences"]
        _chainid = ch_input["protein"]["id"]
        origchainid = pdb.origchainids[_chainid]
        
        keep = pdb.chainkeepindices[origchainid]
        origseq = pdb.chainseqs[origchainid]

        function updatedict!(dict, key, newseq::String)
            newseq = collect(newseq)
            newseq[keep] .= collect(dict[key])
            dict[key] = join(newseq[first(keep):last(keep)])
        end

        updatedict!(ch_input["protein"], "sequence", origseq)
        for i in eachindex(ch_input["protein"]["msa"])
            updatedict!(ch_input["protein"]["msa"], i, String(fill('-', length(origseq))))
        end
    end
    molecularinput
end
fillgaps!(x, pdb::PreprocessedPDB{PDBFormat}; warn_pdbformat=true) = warn_pdbformat && @warn "Will not fill in potential missing gaps before prediction, for PDB format. This is only done with MMCIFFormat"

function resetgaps!(prot::ProteinStructure, pdb::PreprocessedPDB{MMCIFFormat}; kws...)
    for (_chainid, origchainid) in pairs(pdb.origchainids)
        keep = pdb.chainkeepindices[origchainid]
        prefixsize = first(keep) - 1
        keep = keep .- prefixsize

        prot[_chainid] = prot[_chainid][keep]
        prot[_chainid].numbering .= 1:length(prot[_chainid])
    end
end
resetgaps!(x, pdb::PreprocessedPDB{PDBFormat}; warn_pdbformat=true) = warn_pdbformat && @warn "Will not fill in potential missing gaps before prediction, for PDB format. This is only done with MMCIFFormat"

function resetchainids!(prot::ProteinStructure, processed_pdb::PreprocessedPDB)
    for i in eachindex(prot.chains)
        ch = prot.chains[i]
        prot.chains[i] = ProteinChain(processed_pdb.origchainids[ch.id], ch.atoms, ch.sequence)
    end
    sort!(prot.chains, by=ch->ch.id)
end

function resetprocessedpdb(processed_pdb::PreprocessedPDB)
    prot = readpdb(processed_pdb.pdb_filename)
    resetchainids!(prot, processed_pdb)
    for ch in prot.chains
        ch.numbering .= 1:length(ch)
    end
    prot
end

function preprocess(proteinfile::String, dir::String)
    mkpath(dir)
    
    origstruc = readfile(proteinfile)
    origchainids = convertchainids!(origstruc)

    processed_pdb_filename = joinpath(dir, origstruc.name * ".pdb")
    writepdb(processed_pdb_filename, origstruc, backboneselector; expand_disordered=false)

    if fileformat(proteinfile) == MMCIFFormat
        chainseqs, chainkeepindices = find_gapinfo(proteinfile, chainresnumbers(readfile(processed_pdb_filename), origchainids), fileformat(proteinfile))
    elseif fileformat(proteinfile) == PDBFormat
        chainseqs, chainkeepindices = Dict{String, String}(), Dict{String, Vector{Int}}()
    end

    return PreprocessedPDB(processed_pdb_filename, origchainids, chainseqs, chainkeepindices, fileformat(proteinfile))
end

end