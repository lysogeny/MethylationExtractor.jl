"""Genome

A wrapper around the FASTX.FASTA reader that provides memoization.
Provide a filename or fasta reader to construct this. A fasta index is necessary as behaviour without one is undefined.
"""
struct Genome
    reader::FASTX.FASTA.Reader
end
Genome(s::String) = isfile("$s.fai") ? Genome(open(FASTX.FASTA.Reader, s, index = "$s.fai")) : error("Could not locate fasta index at `$s.fai`")

"""read(g::Genome, id::AbstractString)

Fetch sequence with name `id` from genome `g`.
"""
Memoize.@memoize function Base.read(g::Genome, id::AbstractString)
    FASTX.FASTA.sequence(g.reader[id])
end
function Base.read(g::Genome, id::Integer)
    read(g, "$id")
end
