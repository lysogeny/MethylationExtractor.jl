using MethylationExtractor

MethylationExtractor.parse_cli(false)
filter_motif("cCwgg", "test/GRCm39.chr19.fa", "test/full.cov.gz", tempname())
filter_motif("wCg", "test/GRCm39.chr19.fa", "test/full.cov.gz", tempname())
filter_motif("gCh", "test/GRCm39.chr19.fa", "test/full.cov.gz", tempname())
filter_motif("hCh", "test/GRCm39.chr19.fa", "test/full.cov.gz", tempname())
