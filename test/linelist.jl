using EchelleCCFs
using Test

    @testset "Linelists/mask I/O" begin
        mask_fn = joinpath(pkgdir(EchelleCCFs),"data","masks","G2.espresso.mas")
        @test isfile(mask_fn)
        @test_nowarn EchelleCCFs.read_linelist_espresso(mask_fn)
        mask_fn = joinpath(pkgdir(EchelleCCFs),"data","line_lists","VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
        @test isfile(mask_fn)
        @test_nowarn EchelleCCFs.read_mask_vald(mask_fn)
        mask_fn = joinpath(pkgdir(EchelleCCFs),"data","masks","HD101501_f75s75v3_mask.csv")
        @test isfile(mask_fn)
        @test_nowarn EchelleCCFs.read_linelist_rvspectml(mask_fn)
    end
