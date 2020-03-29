# dense Nauty graphs

WORDSIZE == 64 || error("WORDSIZE is not 64. Probably all hell will break loose.")

export DenseNautyGraph, DenseNautyDiGraph

"""
A dense Nauty graph is a bit matrix: if the graph has n vertices,
then it is coded as an (WORDSIZE*m)×n bit matrix, with m=(n+WORDSIZE-1)÷WORDSIZE.
Each column of that matrix is thus a chunk of memory of size m, and the ith column represents the adjacency list of the the ith vertex. The bits are reversed in each word, so to know if the ith vertex is connected to vertex j, check data[setpos(j),i].
"""
struct DenseNautyGraph <: AbstractGraph{Int}
    data::BitMatrix
end
LightGraphs.is_directed(::Type{DenseNautyGraph}) = false
LightGraphs.is_directed(::DenseNautyGraph) = false

struct DenseNautyDiGraph <: AbstractGraph{Int}
    data::BitMatrix
end
LightGraphs.is_directed(::Type{DenseNautyDiGraph}) = true
LightGraphs.is_directed(::DenseNautyDiGraph) = true

DenseNautyXGraph = Union{DenseNautyGraph,DenseNautyDiGraph}

num_setwords(n::Int) = (n+WORDSIZE-1)÷WORDSIZE
setpos(n::Int) = 1 + (((n-1) & ~0x0003f) | (~(n-1) & 0x3f))

DenseNautyGraph(n::Int) = (data = BitMatrix(undef,WORDSIZE*num_setwords(n),n); fill!(data,false); DenseNautyGraph(data))

function DenseNautyGraph(mat::Union{Matrix,BitMatrix})
    n = size(mat,1)
    mat == mat' || error("mat should be a square, symmetric matrix, $mat, $(mat')")

    data = BitMatrix(undef,WORDSIZE*num_setwords(n),n)
    fill!(data,false)
    for s=1:n, d=1:n
        data[setpos(d),s] = mat[s,d]
    end
    DenseNautyGraph(data)
end

function DenseNautyDiGraph(mat::Union{Matrix{Int},Matrix{Bool},BitMatrix})
    n = size(mat,1)
    n == size(mat,2) || error("mat should be a square matrix")

    data = BitMatrix(undef,WORDSIZE*num_setwords(n),n)
    fill!(data,false)
    for s=1:n, d=1:n
        data[setpos(d),s] = mat[s,d]
    end
    DenseNautyDiGraph(data)
end

function DenseNautyGraph(g::AbstractGraph{T}) where T <: Integer
    ng = DenseNautyGraph(nv(g))
    for e=edges(g)
        add_edge!(ng,e)
    end
    ng
end

function DenseNautyDiGraph(g::AbstractGraph{T}) where T <: Integer
    ng = DenseNautyDiGraph(nv(g))
    for e=edges(g)
        add_edge!(ng,e)
    end
    ng
end
        
convert(::Type{DenseNautyGraph}, g) = DenseNautyGraph(g)
convert(::Type{DenseNautyDiGraph}, g) = DenseNautyDiGraph(g)

import Base.==
==(g::DenseNautyXGraph, h::DenseNautyXGraph) = g.data == h.data

Base.eltype(g::DenseNautyXGraph) = Int
LightGraphs.edgetype(g::DenseNautyXGraph) = SimpleEdge{eltype(g)}

LightGraphs.nv(g::DenseNautyXGraph) = size(g.data,2)

LightGraphs.ne(g::DenseNautyGraph) = count(g.data) ÷ 2
LightGraphs.ne(g::DenseNautyDiGraph) = count(g.data)

Base.zero(g::DenseNautyGraph) = DenseNautyGraph(BitMatrix(undef,0,0))
Base.zero(g::DenseNautyDiGraph) = DenseNautyDiGraph(BitMatrix(undef,0,0))

Base.reverse(g::DenseNautyGraph) = g
function Base.reverse(g::DenseNautyDiGraph)
    n = nv(g)
    data = BitMatrix(undef,WORDSIZE*num_setwords(n),n)
    fill!(data,false)
    for s=1:n, d=1:n
        data[setpos(d),s] = g.data[setpos(s),d]
    end
    DenseNautyDiGraph(data)
end

LightGraphs.edges(g::DenseNautyXGraph) = DenseNautyEdgeIter(g)
struct DenseNautyEdgeIter <: AbstractEdgeIter
    g::DenseNautyXGraph
end
function Base.iterate(eit::DenseNautyEdgeIter, state=(1,0))
    g = eit.g
    n = nv(g)
    s, d = state
    while true
        d += 1
        if d > n
            s += 1
            d = 1
        end
        if s > n
            return nothing
        end
        if g.data[setpos(d),s]
            return Edge(s,d), (s,d)
        end
    end
end

LightGraphs.has_edge(g::DenseNautyXGraph, s::Int, d::Int) = g.data[setpos(d),s]

LightGraphs.inneighbors(g::DenseNautyXGraph, v::Int) = g.data[setpos(v),:]

LightGraphs.outneighbors(g::DenseNautyXGraph, v::Int) = g.data[:,v]

LightGraphs.add_edge!(g::DenseNautyGraph, s::Int, d::Int) = (g.data[setpos(d),s] = g.data[setpos(s),d] = true)
LightGraphs.add_edge!(g::DenseNautyDiGraph, s::Int, d::Int) = (g.data[setpos(d),s] = true)
LightGraphs.add_edge!(g::DenseNautyXGraph, e::Edge) = add_edge!(g, e.src, e.dst)

LightGraphs.rem_edge!(g::DenseNautyGraph, s::Int, d::Int) = (g.data[setpos(d),s] = g.data[setpos(s),d] = false; true)
LightGraphs.rem_edge!(g::DenseNautyDiGraph, s::Int, d::Int) = (g.data[setpos(d),s] = false; true)
LightGraphs.rem_edge!(g::DenseNautyXGraph, e::Edge) = rem_edge!(g, e.src, e.dst)
