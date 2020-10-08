using EchelleCCFs
using Test

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
