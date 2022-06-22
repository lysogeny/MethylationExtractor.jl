# MethylationExtractor.jl

A flexible replacement for bismark's `coverage2cytosine` tool.

## Arbitrary context

This tool allows you to define an arbitrary context to filter the `cov` files on.
The context is defined as an IUPAC sequence, the first uppercased letter is understood as the methylated base. 

Examples:

- wCg (CpG methylation in bismark)
- gCh (GpC methylation in bismark)
- cCwgg (DCM methylation)
- asCctws (sequence that has a methylation on the third base)
- ascCtws (same sequence that has a methylation on the fourth base)

# Usage

```
methylation_filter -m MOTIF [-o OUT] fasta in
```

- `fasta` is the reference fasta file. This is required.
- `in` is the input file. This is required.
- `MOTIF` is the motif to be considered. Please see the above section on how to write these.
- `OUT` is an optional output file. If this is not defined the output is printed to `stdout`.

# Installation

Download a copy of the pre-built binaries from the releases page. 

If you want to build this yourself install `PackageCompiler.jl`, get a copy of this repository and use the `scripts/create_image.jl` to create a directory with binaries.

Should you want to use the methods in this package, install this package. In your julia REPL type:

    ] add https://github.com/lysogeny/MethylationExtractor.jl
