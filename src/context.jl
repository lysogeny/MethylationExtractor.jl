"""fetch_context(g::Genome, m::Motifg, ce::CSV.Row)

For the entry `ce` in a coverage file fetch the surrounding bases from genome `g` assuming that the base is part of motif `m`.
Note that this function will return a truncated motif if the position is on the edge of our reference.
"""
function fetch_context(g::Genome, m::Motif, ce::Position)
    start_pos = ce.pos - length_before(m)
    end_pos =  ce.pos + length_after(m)
    sequence = read(g, ce.seq)
    if start_pos < 1
        start_pos = 1
    end
    if end_pos > length(sequence)
        end_pos = length(sequence)
    end
    sequence[start_pos:end_pos]
end

"""iscontext(g::Genome, m::Motifg, ce::CSV.Row)

Return true if entry `ce`'s base fits into the context defined by the motif `m` in the genome `g`.
"""
function iscontext(g::Genome, m::Motif, ce::Position)
    context = fetch_context(g, m, ce)
    for n in findall(context .== BioSequences.DNA_N)
        # Replace the missing bases with gaps for our check.
        # Gaps will not match to our motif.
        # This is the behaviour of our reference implementation (Bismarck)
        context[n] = BioSequences.DNA_Gap
    end
    occursin(BioSequences.ExactSearchQuery(context, BioSequences.iscompatible), m.seq)
end

"""iscontextorreverse(g::Genome, m::Motifg, ce::CSV.Row)

Return true if entry `ce`'s base fits into the context defined by the motif `m` in the genome `g` or the reverse complement.
"""
iscontextorreverse(g::Genome, m::Motif, pos::Position) = iscontext(g, m, pos) | iscontext(g, BioSequences.reverse_complement(m), pos)
iscontextorreverse(g::Genome, m::Motif, ce::CSV.Row) = iscontextorreverse(g, m, Position(ce))
