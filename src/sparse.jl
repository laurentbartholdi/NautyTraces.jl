# sparse Nauty graphs

const sg_weight = Cint

struct SparseNautyGraph <: AbstractGraph{Int}
    nv::Cint # number of vertices
    nde::Cint # number of directed edges (1 for loops)
    v::Ptr{Cint} # for each vertex, its neighbour list is at e[v[i]:v[i]+d[i]-1]
    d::Ptr{Cint} # (out)degrees of vertices
    e::Ptr{Cint} # list of neighourhoods, see above
    w::Ptr{sg_weight} # unused, should be NULL
    vlen::Int # allocated length of v
    dlen::Int # ... d
    elen::Int # ... e
    wlen::Int # ... w; should be 0
end
