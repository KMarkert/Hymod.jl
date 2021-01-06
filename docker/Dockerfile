FROM plutojl/pluto:latest

# Add additional packages
RUN julia -e 'using Pkg; Pkg.add("HTTP"); Pkg.add("CSV"); Pkg.add("DataFrames"); Pkg.add("Plots"); Pkg.add(PackageSpec(url="https://github.com/JuliaHydro/Hymod.jl"))'
