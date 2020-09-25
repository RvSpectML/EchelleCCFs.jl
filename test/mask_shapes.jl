using EchelleCCFs
using Test

    @testset "CCF shape constructors" begin
        Δv = 7e3
        c = EchelleCCFs.speed_of_light_mps
        Δz = Δv/c
        λ = 5000.0

        import EchelleCCFs: λ_min, λ_max, integrate
        @testset "Tophat mask" begin
            @test_nowarn TopHatCCFMask(Δv)
            m = TopHatCCFMask(Δv)
            @test_nowarn λ_min(m,λ)
            @test_nowarn λ_max(m,λ)
            @test_nowarn integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ)
            @test integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ) ≈ 1    rtol = 1e-3
        end
        @testset "Half Cos mask" begin
            @test_nowarn CosCCFMask(Δv)
            m = CosCCFMask(Δv)
            @test_nowarn λ_min(m,λ)
            @test_nowarn λ_max(m,λ)
            @test_nowarn integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ)
            @test integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ) ≈ 1    rtol = 1e-3
        end
        @testset "Gaussian mask" begin
            @test_nowarn GaussianCCFMask(Δv)
            m = GaussianCCFMask(Δv)
            @test_nowarn λ_min(m,λ)
            @test_nowarn λ_max(m,λ)
            @test_nowarn integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ)
            @test integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ) ≈ 1    rtol = 1e-3
        end
        @testset "SuperGaussian mask" begin
            @test_nowarn SuperGaussianCCFMask(Δv,1.1)
            m = SuperGaussianCCFMask(Δv,1.1)
            @test_nowarn λ_min(m,λ)
            @test_nowarn λ_max(m,λ)
            @test_nowarn integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ)
            @test integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ) ≈ 1    rtol = 1e-3
        end
        #=  TODO: Test before using
        @testset "Gaussian Mixture mask" begin
            @test_nowarn GaussianMixtureCCFMask([0.75, 0.25], [Δv, 0.5*Δv], 2*Δv)
            m = GaussianMixtureCCFMask([0.75, 0.25], [Δv, 0.5*Δv], 2*Δv)
            @test_nowarn λ_min(m,λ)
            @test_nowarn λ_max(m,λ)
            @test_nowarn integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ)
            @test integrate(m, (λ_min(m,λ)-λ)*c/λ,  (λ_max(m,λ)-λ)*c/λ) ≈ 1    rtol = 1e-3
        end
        =#
    end
