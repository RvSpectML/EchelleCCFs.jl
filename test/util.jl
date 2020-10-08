using EchelleCCFs
using RvSpectMLBase
using Test

@testset "Utilities" begin
    @testset "Physics" begin
        @test EchelleCCFs.λ_vac_to_air(EchelleCCFs.λ_air_to_vac(5000.0)) ≈ 5000.0
    end

    @testset "Make mask utils" begin
        Δv = 7e3
        Δz = Δv/EchelleCCFs.speed_of_light_mps
        @test_nowarn EchelleCCFs.ChunkWidthFixedΔlnλ(Δz)
    end

end
