using EchelleCCFs
using RvSpectMLBase
using Test

@testset "EchelleCCFs.jl" begin
    @testset "Physics" begin
        @test EchelleCCFs.λ_vac_to_air(EchelleCCFs.λ_air_to_vac(5000.0)) ≈ 5000.0
    end
    @testset "Linelists/mask I/O" begin
        mask_fn = joinpath(pkgdir(EchelleCCFs),"data","masks","G2.espresso.mas")
        @test isfile(mask_fn)
        @test_nowarn EchelleCCFs.read_linelist_espresso(mask_fn)
        @test_nowarn EchelleCCFs.read_mask_espresso(mask_fn)
        # TODO: Add test for VALD once add example data file to repo
    end

    @testset "Make mask utils" begin
        Δv = 7e3
        Δz = Δv/EchelleCCFs.speed_of_light_mps
        @test_nowarn EchelleCCFs.ChunkWidthFixedΔlnλ(Δz)
    end

    include("mask_shapes.jl")

    @testset "CCF plan" begin
        Δv = 7e3
        λ1 = 5025.0
        d1 = 0.5
        d2 = 0.25
        λ2 = 5075.0
        @test_nowarn BasicLineList([λ1,λ2], [d1,d2])
        ll = BasicLineList([λ1,λ2], [d1,d2])
        @test_nowarn BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        ccfpl = BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        @test_nowarn calc_ccf_v_grid(ccfpl)
        v_grid = calc_ccf_v_grid(ccfpl)
        @test length(v_grid) == calc_length_ccf_v_grid(ccfpl)
    end

    include("calc_rv.jl")

    include("convenience.jl")

end
