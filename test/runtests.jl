using Graphs
using NautyTraces
using Test

@testset "NautyTraces.jl" begin
    # Petersen graph
    g = DenseNautyGraph(10)
    for i=1:5
        add_edge!(g,i,mod1(i+1,5))
        add_edge!(g,i,5+i)
        add_edge!(g,5+i,5+mod1(i+2,5))
    end
    @test nv(g) == 10
    @test ne(g) == 15
    result = nauty(g,automgroup=true)
    @test result.grpsize == 120
end
