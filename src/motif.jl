"""Motif

encodes information about the methylation motif.
The first uppercase character is understood as the methylation site.
motif"gC" is a motif where a G precedes a methylated C.

Examples:

    julia> motif"cCwgg"
    Motif(CCWGG, 2)

    julia> Motif("cCwgg")
    Motif(CCWGG, 2)

    julia> motif"Cg"
    Motif(CG, 1)

    julia> motif"gC"
    Motif(gC, 2)

"""
struct Motif
    seq::BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}
    pos::Int16
end
function Motif(s::AbstractString)
    chars = collect.(s)
    pos = findfirst(isuppercase.(chars))
    Motif(BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}(s), pos)
end
macro motif_str(s) 
    Motif(s)
end

# Other library extensions.
function Base.length(m::Motif)
    length(m.seq)
end

function BioSequences.reverse_complement(m::Motif)
    Motif(BioSequences.reverse_complement(m.seq), length(m) - length_before(m))
end

"""length_before(m::Motif)

Bases before methylation site
"""
function length_before(m::Motif)
    m.pos-1
end

"""length_after(m::Motif)

Bases after methylation site
"""
function length_after(m::Motif)
    length(m) - m.pos
end

