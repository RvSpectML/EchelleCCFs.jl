using EchelleCCFs
using Test
import Statistics: mean

@testset "Check CCF accuracy" begin  # TODO repeat for multiple mask shapes
    @testset "Tophat mask w/ 1 line" begin
        Δv = 1000
        λ1 = 5025.0
        d1 = 0.5
        σ_v = 5000.0
        λlo = 5000
        λhi = 5050
        npix = 2000
        @test_nowarn BasicLineList([λ1], [d1])
        ll = BasicLineList([λ1], [d1])
        ll = BasicLineList([λ1], [d1])
        @test_nowarn BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        ccfpl = BasicCCFPlan(line_list=ll, mask_shape=TopHatCCFMask(Δv))
        @test_nowarn calc_ccf_v_grid(ccfpl)
        v_grid = calc_ccf_v_grid(ccfpl)
        @test length(v_grid) == calc_length_ccf_v_grid(ccfpl)
        λs = exp.(range(log(λlo),stop=log(λhi),length=npix))
        flux = ones(size(λs))
        flux .*= 1 .- d1.*exp.(-0.5.*((λs.-λ1)./λ1.*(EchelleCCFs.speed_of_light_mps./σ_v)).^2)
        @test_nowarn ccf_1D(λs, flux, ccfpl)
        ccf = ccf_1D(λs, flux, ccfpl)
        @test sum(ccf) > 0
        @test_nowarn  MeasureRvFromCCFGaussian()
        mrfcg = MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.75)
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=mrfcg )
        vfit = measure_rv_from_ccf(v_grid,ccf, alg=mrfcg )
        @test vfit.rv ≈ 0  atol = 0.1   # TODO: tune tolerance better
        var = flux./10^5
        @test_nowarn ccf_1D(λs, flux, var, ccfpl)
        (ccf, ccf_var) = ccf_1D(λs, flux, var, ccfpl)
        @test sum(ccf) > 0
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, ccf_var, alg=mrfcg )
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=mrfcg )  # TODO: Check Why is this doing so much worse?
        #=
        using Plots
        global plt = plot()
        plot!(plt,v_grid,ccf,label=:none)
        scatter!(plt,v_grid,ccf,color=1,label=:none)
        plot!(plt,v_grid,sqrt.(ccf_var),label=:none)
        =#
        @test vfit.rv ≈ 0  atol = 0.5   # TODO: tune tolerance better
    end

    @testset "Tophat mask w/ 2 lines" begin
        Δv = 1000
        λ1 = 5025.0
        λ2 = 5075.0
        d1 = 0.5
        d2 = 0.25
        σ_v = 5000.0
        λlo = 5000
        λhi = 5100
        npix = 8000
        @test_nowarn BasicLineList([λ1,λ2], [d1,d2])
        ll = BasicLineList([λ1,λ2], [d1,d2])
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
        @test vfit.rv ≈ 0  atol = 0.3   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=mrfcg )  # TODO: Check Why is this doing so much worse?
        #=
        using Plots
        global plt = plot()
        plot!(plt,v_grid,ccf,label=:none)
        scatter!(plt,v_grid,ccf,color=1,label=:none)
        plot!(plt,v_grid,sqrt.(ccf_var),label=:none)
        =#
        @test vfit.rv ≈ 0  atol = 0.3   # TODO: tune tolerance better

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
        npix = 4000
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
        @test vfit.rv ≈ 0  atol = 0.3   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=MeasureRvFromCCFQuadratic() )
        @test vfit.rv ≈ 0  atol = 0.3   # TODO: tune tolerance better

        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromMinCCF() )
        vfit = measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromMinCCF())
        @test vfit.rv ≈ 0  atol = 0.1   # TODO: tune tolerance better
        vfit = measure_rv_from_ccf(v_grid,ccf,ccf_var, alg=MeasureRvFromMinCCF())
        @test vfit.rv ≈ 0  atol = 0.5   # TODO: tune tolerance better

        #=
        @test_nowarn measure_rv_from_ccf(v_grid,ccf, alg=MeasureRvFromCCFCentroid() )
        vfit = measure_rv_from_ccf(v_grid,ccf)
        @test vfit.rv ≈ 0  atol = 500   # TODO: tune tolerance better
        =#
    end
    @testset "RV injection/extraction" begin
        # velocities
        amp = 1.0 # meters per second
        vel = amp .* sin.(range(0, 2π, length=100))
        mid = 5434.5
        dep = 0.75
        wid = 0.1
        mid_shift = mid .* EchelleCCFs.calc_doppler_factor.(vel)

        # compute gaussian line spectrum
        buffer = 0.75
        resolution = 7e5
        Δlnλ = (1.0 / resolution)
        wavs = exp.(range(log(mid - buffer), log(mid + buffer), step=Δlnλ))

        flux = zeros(length(wavs), length(vel))
        for i in eachindex(vel)
            flux[:,i] .= @. 1.0 - dep * exp(-((wavs-mid_shift[i])/wid)^2 / 2.0)
        end

        # measure ccf
        line_list = EchelleCCFs.BasicLineList([mid], [dep])
        speed_of_light = EchelleCCFs.speed_of_light_mps
        mask_width = speed_of_light/resolution
        mask_shape = TopHatCCFMask(mask_width)
        ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape)
        v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
        ccf = EchelleCCFs.ccf_1D(wavs, flux, ccf_plan)

        # measure velocities from CCF
        mrv_gauss = MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.5)
        mrv_quad = MeasureRvFromCCFQuadratic(frac_of_width_to_fit=0.1)

        rvs_gauss = zeros(length(vel))
        sig_gauss = zeros(length(vel))
        rvs_quad = zeros(length(vel))
        sig_quad = zeros(length(vel))
        for i in eachindex(vel)
            rvs_gauss[i], sig_gauss[i] = mrv_gauss(v_grid, view(ccf, :, i))
            rvs_quad[i], sig_quad[i] = mrv_quad(v_grid, view(ccf, :, i))
        end
        err_gauss = abs.(rvs_gauss .- mean(rvs_gauss) .- vel)
        err_quad = abs.(rvs_quad .- mean(rvs_quad) .- vel)

        # test maximum errors are less than a cm/s
        @assert maximum(err_gauss) < 1e-2
        @assert maximum(err_quad) < 1e-2
    end
end
