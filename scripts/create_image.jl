using PackageCompiler
using Pkg
Pkg.activate(".")
create_app(".", "build", executables=["methylation_filter" => "julia_main"], force=true, precompile_statements_file="precompile_statements.jl")
