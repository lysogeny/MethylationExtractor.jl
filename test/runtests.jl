using BioSequences
using MethylationExtractor
using CSV

using Test

@testset "Motif" begin
    @test length(motif"wCg") == 3
    @test length_before(motif"wCg") == 1
    @test length_after(motif"wCg") == 1
    @test length_after(motif"wCgggg") == 4
    @test length_after(motif"wwwwCgggg") == 4
    @test length(motif"wwwwCgggg") == 9
    m1 = motif"wCg"
    m2 =  Motif("wCg")
    @test m1.seq == m2.seq
    @test m1.pos == m2.pos
end

@testset "Genome" begin
    genome = Genome("GRCm39.chr19.fa")
    seq1 = read(genome, "19")
    @test length(seq1) > 0
    seq2 = read(genome, 19)
    @test length(seq2) > 0
    seq1 == seq2
end

@testset "context fetching" begin
    motif = motif"wCg"
    genome = Genome("GRCm39.chr19.fa")
    cpg_file = "NOMe.CpG.cov.gz"
    #dcm_file = "DCM.cov.gz"
    cpg_entries = CSV.File(cpg_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    # Test all of them contained in the correct thing
    #println(fetch_context(genome, motif_cpg, cpg_entries[1]))
    @test iscontextorreverse(genome, motif, cpg_entries[1])
    #@test all(map(x -> iscontextorreverse(genome, motif_cpg, x), cpg_entries))
end

@testset "sample files in context" begin
    motif_cpg = motif"wCg"
    motif_gpc = motif"gCh"
    genome = Genome("GRCm39.chr19.fa")
    cpg_file = "NOMe.CpG.cov.gz"
    gpc_file = "NOMe.GpC.cov.gz"
    cpg_entries = CSV.File(cpg_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    gpc_entries = CSV.File(gpc_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    # Test all of them contained in the correct file
    @test all(map(x -> iscontextorreverse(genome, motif_cpg, x), cpg_entries))
    @test all(map(x -> iscontextorreverse(genome, motif_gpc, x), gpc_entries))
    # Crosstest, CpG motifs are not GpC motifs and GpC motifs are not CpG motifs.
    @test all(map(x -> !iscontextorreverse(genome, motif_gpc, x), cpg_entries))
    @test all(map(x -> !iscontextorreverse(genome, motif_cpg, x), gpc_entries))
end

@testset "Tool internal function works for CpG" begin
    in_file = "full.cov.gz"
    ref_file = "NOMe.CpG.cov.gz"
    out_file = tempname()
    genome_file = "GRCm39.chr19.fa"
    motif = "wCg"
    filter_motif(motif, genome_file, in_file, out_file)
    # Check result makes sense
    result = CSV.File(out_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    ref = CSV.File(ref_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    result_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in result])
    ref_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in ref])
    @test length(result_pos) == length(ref_pos)
    @test length(result_pos) == length(intersect(result_pos, ref_pos))
    @test length(setdiff(result_pos, ref_pos)) == 0
    @test length(setdiff(ref_pos, result_pos)) == 0
end

@testset "Tool internal function works for GpC" begin
    in_file = "full.cov.gz"
    ref_file = "NOMe.GpC.cov.gz"
    out_file = tempname()
    genome_file = "GRCm39.chr19.fa"
    motif = "gCh"
    filter_motif(motif, genome_file, in_file, out_file)
    # Check result makes sense
    result = CSV.File(out_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    ref = CSV.File(ref_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    result_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in result])
    ref_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in ref])
    @test length(result_pos) == length(ref_pos)
    @test length(result_pos) == length(intersect(result_pos, ref_pos))
    @test length(setdiff(result_pos, ref_pos)) == 0
    @test length(setdiff(ref_pos, result_pos)) == 0
end

@testset "Tool internal function works for DCM" begin
    in_file = "full.cov.gz"
    ref_file = "DCM.cov.gz"
    out_file = tempname()
    genome_file = "GRCm39.chr19.fa"
    motif = "cCwgg"
    filter_motif(motif, genome_file, in_file, out_file)
    # Check result makes sense
    result = CSV.File(out_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    ref = CSV.File(ref_file, header=[:seq, :pos_start, :pos_stop, :perc_meth, :n_meth, :n_unmeth])
    result_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in result])
    ref_pos = Set(["$(x.seq) $(x.pos_start):$(x.pos_stop)" for x in ref])
    @test length(result_pos) == length(ref_pos)
    @test length(result_pos) == length(intersect(result_pos, ref_pos))
    @test length(setdiff(result_pos, ref_pos)) == 0
    @test length(setdiff(ref_pos, result_pos)) == 0
end
