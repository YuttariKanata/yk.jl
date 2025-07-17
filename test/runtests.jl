using yk
using Test

@testset "yk.jl" begin
    @test iroot(1331, 3) == 11
    @test iroot(1331, 4) == 6
    @test iroot(BigInt(2)^1000, 1000) == 2
    @test_throws DomainError iroot(-4, 2)
end
