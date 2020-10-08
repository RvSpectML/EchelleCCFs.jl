using EchelleCCFs
using RvSpectMLBase
using Test

@testset "EchelleCCFs.jl" begin

    include("util.jl")

    include("linelist.jl")

    include("mask_shapes.jl")

    include("ccf_plan.jl")

    include("calc_rv.jl")

    include("convenience.jl")

end
