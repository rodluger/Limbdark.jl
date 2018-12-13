ENV["PYTHON"]=""
using Pkg
Pkg.add("Conda")
using Conda
Conda.update()
Conda.add("matplotlib")
Conda.add_channel("conda-forge")
Conda.add("starry")
Pkg.add("PyCall")
Pkg.add("Cubature")
Pkg.build("PyCall")
Pkg.add("PyPlot")
Pkg.add("SpecialFunctions")
Pkg.add("ForwardDiff")
Pkg.add("DiffResults")
Pkg.add("Optim")
Pkg.add("QuadGK")
