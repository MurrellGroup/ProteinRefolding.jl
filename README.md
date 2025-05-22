# ProteinRefolding

[![Build Status](https://github.com/MurrellGroup/ProteinRefolding.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/ProteinRefolding.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Quick start

```julia
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/MurrellGroup/ProteinRefolding.jl"))

using ProteinRefolding

# Single sequence refolding using ProteinMPNN and Boltz1x
refold("myprotein.cif")
# MSA-refolding using ProteinMPNN and Boltz1x
refold("myprotein.cif", msa_size=100)
```