using EchelleCCFs
using Test


@testset "Convenience functions" begin  # TODO repeat for multiple mask shapes
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
    using RvSpectMLBase
    @testset "calc_ccf_chunk" begin
        chunk = ChunkOfSpectrum(λs,flux,ones(size(flux)))
        @test_nowarn calc_ccf_chunk(chunk, ccfpl)
    end
    @testset "calc_ccf_chunklist" begin
        chunk = ChunkOfSpectrum(λs,flux,ones(size(flux)))
        chunklist = ChunkList([chunk,chunk,chunk],[1,2,3])
        @test_nowarn calc_ccf_chunklist(chunklist, fill(ccfpl,3) )
    end
    @testset "calc_order_ccfs_chunklist" begin
        chunk = ChunkOfSpectrum(λs,flux,ones(size(flux)))
        chunklist = ChunkList([chunk,chunk,chunk],[1,2,3])
        @test_nowarn EchelleCCFs.calc_order_ccfs_chunklist(chunklist, fill(ccfpl,3)  )
    end

    #= Not working
    @testset "calc_ccf_chunklist_timeseries" begin
        chunk = ChunkOfSpectrum(λs,flux,ones(size(flux)))
        chunklist = ChunkList([chunk,chunk,chunk],[1,2,3])
        #clt = ChunkListTimeseries(rand(3), [chunklist,chunklist,chunklist] )
        #@test_nowarn calc_ccf_chunklist_timeseries(clt, fill(ccfpl,3)  )
    end
    =#

    #=  Not working
    @testset "calc_ccf_chunklist_timeseries" begin
        using RvSpectMLBase.TheoreticalInstrument
        using DataFrames
        num_pixels = 1000
        num_orders = 2
        num_obs = 3
        lambda2d = collect(reshape(range(λlo,stop=λhi,length=num_pixels*num_orders),num_pixels,num_orders))
        inst = TheoreticalInstrument2D(lambda2d)
        times = sort(rand(num_obs))
        rvs = range(0.0, stop=10.0, length=length(times))
        line_list = DataFrame(:lambda=>[λ1,λ2], :weight=>[d1,d2])
        @test_nowarn TheoreticalInstrument.generate_spectra_timeseries(times,line_list,inst, rvs)
        spectra = TheoreticalInstrument.generate_spectra_timeseries(times,line_list,inst, rvs)
        #chunk_list = DataFrame(:lambda_lo=>[λlo, (λlo+λhi)/2], :lambda_hi=>[(λlo+λhi)/2, λhi])
        chunk_list = DataFrame(:lambda_lo=>0.9.*[λ1, λ2], :lambda_hi=>1.1.*[λ1, λ2])
        @test_nowarn make_chunk_list_timeseries(spectra,chunk_list)
        clt = make_chunk_list_timeseries(spectra,chunk_list)
        @test_nowarn calc_ccf_chunklist_timeseries(clt, fill(ccfpl,num_obs)  )
    end
    =#

end
