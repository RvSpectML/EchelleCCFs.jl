using EchelleCCFs
using Test


import Polynomials
function calc_vertex_of_quadatic_fit(x::AbstractVector{T1}, y::AbstractVector{T2}) where {T1<:Real, T2<:Real}
    pfit = Polynomials.fit(x, y, 2)
    @assert length(Polynomials.coeffs(pfit)) >= 3   # just in case fails to fit a quadratic
    c, b, a = Polynomials.coeffs(pfit)
    v_at_min_of_quadratic = -b/(2*a)
    return v_at_min_of_quadratic
end

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
        ccf = ccf_1D(λs, flux, ccfpl)
        idx_min = findmin(ccf)[2]
        idx_near_min = idx_min-3:idx_min+3
        vfit = calc_vertex_of_quadatic_fit(v_grid[idx_near_min], ccf[idx_near_min])
        #println("TopHatCCFMask:  vfit = ", vfit)
        @test vfit ≈ 0  atol = 500   # TODO: tune tolerance better
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
        ccf = ccf_1D(λs, flux, ccfpl)
        idx_min = findmin(ccf)[2]
        idx_near_min = idx_min-3:idx_min+3
        vfit = calc_vertex_of_quadatic_fit(v_grid[idx_near_min], ccf[idx_near_min])
        #println("GaussianCCFMask:  vfit = ", vfit)
        @test vfit ≈ 0  atol = 100   # TODO: tune tolerance better
    end
end
