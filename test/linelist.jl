using EchelleCCFs
using Test

    @testset "Linelists/mask I/O" begin
        mask_fn = joinpath(pkgdir(EchelleCCFs),"data","masks","G2.espresso.mas")
        @test isfile(mask_fn)
        @test_nowarn EchelleCCFs.read_linelist_espresso(mask_fn)
        @test_nowarn EchelleCCFs.read_mask_espresso(mask_fn)
        # TODO: Add test for VALD once add example data file to repo
    end
