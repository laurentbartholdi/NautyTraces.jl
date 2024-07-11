# dense Nauty graphs

WORDSIZE == 64 || error("WORDSIZE is not 64. Probably all hell will break loose.")

export DenseNautyGraph, DenseNautyDiGraph

"""
  DenseNautyGraph(size::Int)\\
  DenseNautyGraph(adj::Matrix)\\
  DenseNautyGraph(g::AbstractGraph)

Constructs a new dense Nauty graph.

A dense Nauty graph is a bit matrix: if the graph has `n` vertices,
then it is coded as an `(WORDSIZE*m)×n` bit matrix, with `m=(n+WORDSIZE-1)÷WORDSIZE`.
Each column of that matrix is thus a chunk of memory of size `m`, and the `i`th column represents the adjacency list of the the `i`th vertex. The bits are reversed in each word, so to know if the `i`th vertex is connected to vertex `j`, check the field `data[setpos(j),i]`.
"""
struct DenseNautyGraph <: AbstractGraph{Int}
    data::BitMatrix
end
Graphs.is_directed(::Type{DenseNautyGraph}) = false
Graphs.is_directed(::DenseNautyGraph) = false

"""
  DenseNautyDiGraph(size::Int)\\
  DenseNautyDiGraph(adj::Matrix)\\
  DenseNautyDiGraph(g::AbstractGraph)

Constructs a new dense Nauty directed graph.

A dense Nauty directed graph is a bit matrix: if the graph has `n` vertices,
then it is coded as an `(WORDSIZE*m)×n` bit matrix, with `m=(n+WORDSIZE-1)÷WORDSIZE`.
Each column of that matrix is thus a chunk of memory of size `m`, and the `i`th column represents the adjacency list of the the `i`th vertex. The bits are reversed in each word, so to know if the `i`th vertex is connected to vertex `j`, check the field `data[setpos(j),i]`.
"""
struct DenseNautyDiGraph <: AbstractGraph{Int}
    data::BitMatrix
end
Graphs.is_directed(::Type{DenseNautyDiGraph}) = true
Graphs.is_directed(::DenseNautyDiGraph) = true

DenseNautyXGraph = Union{DenseNautyGraph,DenseNautyDiGraph}

num_setwords(n) = (n+WORDSIZE-1)÷WORDSIZE
setpos(n) = 1+((n-1)⊻0x3f)

DenseNautyGraph(n::T) where T <: Integer = (data = BitMatrix(undef,WORDSIZE*num_setwords(n),n); fill!(data,false); DenseNautyGraph(data))

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

DenseNautyDiGraph(n::T) where T <: Integer = (data = BitMatrix(undef,WORDSIZE*num_setwords(n),n); fill!(data,false); DenseNautyDiGraph(data))

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

Base.:(==)(g::DenseNautyXGraph, h::DenseNautyXGraph) = g.data == h.data
Base.hash(g::DenseNautyXGraph) = hash(g.data)
Base.copy(g::DenseNautyXGraph) = typeof(g)(copy(g.data))
    
Base.eltype(g::DenseNautyXGraph) = Int
Graphs.edgetype(g::DenseNautyXGraph) = Graphs.SimpleGraphs.SimpleEdge{eltype(g)}

Graphs.nv(g::DenseNautyXGraph) = size(g.data,2)

function Graphs.ne(g::DenseNautyGraph)
    n = 0
    for i=1:nv(g), j=1:i
        if g.data[setpos(i),j] || g.data[setpos(j),i]
            n += 1
        end
    end
    n
end
Graphs.ne(g::DenseNautyDiGraph) = count(g.data)

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

Graphs.edges(g::DenseNautyGraph) = DenseNautyEdgeIter(g)
struct DenseNautyEdgeIter <: AbstractEdgeIter
    g::DenseNautyGraph
end
Base.eltype(eit::DenseNautyEdgeIter) = Graphs.SimpleGraphs.SimpleEdge{eltype(eit.g)}
Base.length(eit::DenseNautyEdgeIter) = ne(eit.g)
function Base.iterate(eit::DenseNautyEdgeIter, state=(0,1))
    g = eit.g
    n = nv(g)
    d, s = state
    while true
        d += 1
        if d >= s
            s += 1
            d = 1
        end
        if s > n
            return nothing
        end
        if g.data[setpos(d),s]
            return Edge(d,s), (d,s)
        end
    end
end

Graphs.edges(g::DenseNautyDiGraph) = DenseNautyDiEdgeIter(g)
struct DenseNautyDiEdgeIter <: AbstractEdgeIter
    g::DenseNautyDiGraph
end
Base.eltype(eit::DenseNautyDiEdgeIter) = Graphs.SimpleGraphs.SimpleEdge{eltype(eit.g)}
Base.length(eit::DenseNautyDiEdgeIter) = ne(eit.g)
function Base.iterate(eit::DenseNautyDiEdgeIter, state=(1,0))
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

Graphs.has_edge(g::DenseNautyXGraph, s, d) = g.data[setpos(d),s]

function Graphs.inneighbors(g::DenseNautyXGraph, v)
    n = Int[]
    for i=1:nv(g)
        g.data[setpos(v),i] && push!(n,i)
    end
    n
end

function Graphs.outneighbors(g::DenseNautyXGraph, v)
    n = Int[]
    for i=1:nv(g)
        g.data[setpos(i),v] && push!(n,i)
    end
    n
end

Graphs.add_edge!(g::DenseNautyGraph, s, d) = (g.data[setpos(d),s] = g.data[setpos(s),d] = true)
Graphs.add_edge!(g::DenseNautyDiGraph, s, d) = (g.data[setpos(d),s] = true)
Graphs.add_edge!(g::DenseNautyXGraph, e::Edge) = add_edge!(g, e.src, e.dst)

Graphs.rem_edge!(g::DenseNautyGraph, s, d) = (g.data[setpos(d),s] = g.data[setpos(s),d] = false; true)
Graphs.rem_edge!(g::DenseNautyDiGraph, s, d) = (g.data[setpos(d),s] = false; true)
Graphs.rem_edge!(g::DenseNautyXGraph, e::Edge) = rem_edge!(g, e.src, e.dst)
