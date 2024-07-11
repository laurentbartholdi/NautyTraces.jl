# NautyTraces.jl

This package provides Julia bindings to the nauty and traces library, accessible via the nauty graph format or any SimpleGraph format.

```julia
julia> using Pkg; Pkg.add("NautyTraces"); Pkg.build("NautyTraces")

julia> using Graphs, NautyTraces

julia> g = DenseNautyGraph(10) # a graph on 10 vertices
DenseNautyGraph(Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0])

julia> for i=1:5 # the Petersen graph
           add_edge!(g,i,mod1(i+1,5))
           add_edge!(g,i,5+i)
           add_edge!(g,5+i,5+mod1(i+2,5))
       end

julia> result = nauty(g,automgroup=true)
(orbits: IntDisjointSets{Int64}([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), group size: 120, generators: Permutations.Permutation[(1)(2)(3)(4,8)(5,6)(7)(9,10), (1)(2)(3,7)(4,9)(5,6)(8,10), (1)(2,5)(3,4)(6)(7,10)(8,9), (1,2)(3,5)(4)(6,7)(8,10)(9)])

julia> g = DenseNautyGraph(cycle_graph(10)) # a cycle
DenseNautyGraph(Bool[0 0 … 0 0; 0 0 … 0 0; … ; 1 0 … 0 0; 0 1 … 0 1])

julia> result = nauty(g,getcanon=true)
(orbits: DataStructures.IntDisjointSets{Int64}([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), group size: 20, labeling: (1)(2,10,9,3)(4,6,7)(5)(8), canong)

julia> SimpleGraph(result.canong)
{10, 10} undirected simple Int64 graph

julia> collect(edges(ans)) # edges in the canonical ordering of the cycle
10-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 1 => 3
 Edge 2 => 10
 Edge 3 => 9
 Edge 4 => 5
 Edge 4 => 6
 Edge 5 => 7
 Edge 6 => 8
 Edge 7 => 9
 Edge 8 => 10
