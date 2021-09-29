using LightGraphs
using Nauty
using Test

@testset "Nauty.jl" begin
    # Petersen graph
    g = DenseNautyGraph(10)
    for i=1:5
        add_edge!(g,i,1+(i%5))
        add_edge!(g,i,i+5)
        add_edge!(g,i+5,6+((i+2)%5))
    end
    @test nv(g) == 10
    @test ne(g) == 15
    dict = nauty(g,automgroup=true)
    @test dict[:grpsize] == 120
end
