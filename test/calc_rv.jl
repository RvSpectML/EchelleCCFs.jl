using EchelleCCFs
using Test


@testset "Check CCF accuracy" begin  # TODO repeat for multiple mask shapes
    @testset "Tophat mask" begin
        Δv = 1000
        λ1 = 5025.0
        λ2 = 5075.0
        d1 = 0.5
        d2 = 0.25
        σ_v = 5000.0
        λlo = 5000
        λhi = 5100
        npix = 1000
        @test_nowarn BasicLineList([λ1,λ2], [d1,d2])
        ll = BasicLineList([λ1,λ2], [d1,d2])
        @test_nowarn BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        ccfpl = BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        @test_nowarn calc_ccf_v_grid(ccfpl)
        v_grid = calc_ccf_v_grid(ccfpl)
        @test length(v_grid) == calc_length_ccf_v_grid(ccfpl)
        λs = exp.(range(log(λlo),stop=log(λhi),length=npix))
        flux = ones(size(λs))
        flux .*= 1 .- d1.*exp.(-0.5.*((λs.-λ1)./λ1.*(EchelleCCFs.speed_of_light_mps./σ_v)).^2)
        flux .*= 1 .- d2.*exp.(-0.5.*((λs.-λ2)./λ2.*(EchelleCCFs.speed_of_light_mps./σ_v)).^2)
        @test_nowarn ccf_1D(λs, flux, ccfpl)
        var = flux./10^5
        @test_nowarn ccf_1D(λs, flux, var, ccfpl)
        (ccf, ccf_var) = ccf_1D(λs, flux, var, ccfpl)

        @test_nowarn  MeasureRvFromCCFGaussian()
        mrfcg = MeasureRvFromCCFGaussian()
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=mrfcg )
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, ccf_var, alg=mrfcg )
        vfit = measure_rv_from_ccf(v_grid,ccf, alg=mrfcg )
        @test vfit.rv ≈ 0  atol = 0.1   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=mrfcg )  # TODO: Check Why is this doing so much worse?
        @test vfit.rv ≈ 0  atol = 3.0   # TODO: tune tolerance better

    end
    @testset "Gaussian mask" begin
        Δv = 7e3
        λ1 = 5025.0
        λ2 = 5075.0
        d1 = 0.5
        d2 = 0.25
        σ_v = 5000.0
        λlo = 5000
        λhi = 5100
        npix = 1000
        @test_nowarn BasicLineList([λ1,λ2], [d1,d2])
        ll = BasicLineList([λ1,λ2], [d1,d2])
        @test_nowarn BasicCCFPlan(line_list=ll, mask_shape=GaussianCCFMask(σ_v))
        ccfpl = BasicCCFPlan(line_list=ll, mask_shape=GaussianCCFMask(σ_v))
        @test_nowarn calc_ccf_v_grid(ccfpl)
        v_grid = calc_ccf_v_grid(ccfpl)
        @test length(v_grid) == calc_length_ccf_v_grid(ccfpl)
        λs = exp.(range(log(λlo),stop=log(λhi),length=npix))
        flux = ones(size(λs))
        flux .*= 1 .- d1.*exp.(-0.5.*((λs.-λ1)./λ1.*(EchelleCCFs.speed_of_light_mps./σ_v)).^2)
        flux .*= 1 .- d2.*exp.(-0.5.*((λs.-λ2)./λ2.*(EchelleCCFs.speed_of_light_mps./σ_v)).^2)
        @test_nowarn ccf_1D(λs, flux, ccfpl)
        var = flux./10^5
        @test_nowarn ccf_1D(λs, flux, var, ccfpl)
        (ccf, ccf_var) = ccf_1D(λs, flux, var, ccfpl)

        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromCCFQuadratic() )
        @test_nowarn measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=MeasureRvFromCCFQuadratic() )
        vfit = measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromCCFQuadratic() )    # TODO: Check Why is this doing so much worse?
        @test vfit.rv ≈ 0  atol = 3.0   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=MeasureRvFromCCFQuadratic() )
        @test vfit.rv ≈ 0  atol = 3.0   # TODO: tune tolerance better

        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromMinCCF() )
        vfit = measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromMinCCF())
        @test vfit.rv ≈ 0  atol = 0.1   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=MeasureRvFromMinCCF())
        @test vfit.rv ≈ 0  atol = 1.0   # TODO: tune tolerance better

        #=
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromCCFCentroid() )
        vfit = measure_rv_from_ccf(v_grid,ccf)
        @test vfit.rv ≈ 0  atol = 500   # TODO: tune tolerance better
        =#
    end
end
