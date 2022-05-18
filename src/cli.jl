function parse_cli()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "fasta"
            help = "Reference fasta (with accompanying index)"
            required = true
        "in"
            help = "Input coverage file"
            required = true
        "--motif", "-m"
            help = "Sequence motif to filter for"
            required = true
        "--out", "-o"
            help = "Output filtered coverage file"
    end
    ArgParse.parse_args(ARGS, s)
end

function filter_motif(motif::AbstractString, genome::AbstractString, in::AbstractString)
    motif = Motif(motif)
    genome = Genome(genome)
    file = CSV.File(in, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    filter(x -> iscontextorreverse(genome, motif, x), file)
end

function filter_motif(motif::AbstractString, genome::AbstractString, in::AbstractString, out::AbstractString)
    entries = filter_motif(motif, genome, in)
    entries |> CSV.write(out, delim="\t", header=false)
end

function filter_motif(motif::AbstractString, genome::AbstractString, in::AbstractString, out::Nothing)
    entries = filter_motif(motif, genome, in)
    map(print, collect(CSV.RowWriter(entries, delim="\t"))[2:end])
end

function julia_main()::Cint
    args = parse_cli()
    filter_motif(args["motif"], args["fasta"], args["in"], args["out"])
    return 0
end
