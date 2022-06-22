module MethylationExtractor

import FASTX
import BioSequences
import Memoize
import CSV
import ArgParse

include("motif.jl")
include("genome.jl")
include("positions.jl")
include("context.jl")
include("cli.jl")

export Motif, Genome, length_before, length_after, fetch_context, iscontext, iscontextorreverse, @motif_str, filter_motif

end # module
