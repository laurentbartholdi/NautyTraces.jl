module NautyTraces

using Graphs
using Permutations
using Libdl
using DataStructures

export nauty, traces, add_edge!, NautyIsomorphism

const LIB_FILE = "$(@__DIR__)" * "/../deps/nauty." * Libdl.dlext

const WORDSIZE = ccall((:wordsize, LIB_FILE), Int, ())

include("dense.jl")
include("sparse.jl")

# data structure copied from https://github.com/bovine3dom/Nauty.jl

const Nbool = Cint
mutable struct optionblk
    getcanon::Cint	# make canong and canonlab?
    digraph::Nbool	# multiple edges or loops?
    writeautoms::Nbool	# write automorphisms?
    writemarkers::Nbool	# write stats on pts fixed, etc.?
    defaultptn::Nbool	# set lab,ptn,active for single cell?
    cartesian::Nbool	# use cartesian rep for writing automs?
    linelength::Cint	# max chars/line (excl. '\n') for output
    outfile::Ptr{Cvoid}	# FILE *outfile;
		# file for output, if any
    userrefproc::Ptr{Cvoid}	# void (*userrefproc)
		# replacement for usual refine procedure
		# (graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    userautomproc::Ptr{Cvoid}	# void (*userautomproc)
		# procedure called for each automorphism
		# (int,int*,int*,int,int,int);
    userlevelproc::Ptr{Cvoid}	# void (*userlevelproc)
		# procedure called for each level
		# (int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
    usernodeproc::Ptr{Cvoid}	# void (*usernodeproc)
		# procedure called for each node
		# (graph*,int*,int*,int,int,int,int,int,int);
    usercanonproc::Ptr{Cvoid}	# Cint  (*usercanonproc)
		# procedure called for better labellings
		# (graph*,int*,graph*,int,int,int,int);
    invarproc::Ptr{Cvoid}	# void (*invarproc)
		# procedure to compute vertex-invariant
		# (graph*,int*,int*,int,int,int,int*,int,Nbool,int,int);
    tc_level::Cint	# max level for smart target cell choosing
    mininvarlevel::Cint	# min level for invariant computation
    maxinvarlevel::Cint	# max level for invariant computation
    invararg::Cint	# value passed to (*invarproc)()
    dispatch::Ptr{Cvoid}	# dispatchvec *dispatch;
		# vector of object-specific routines
    schreier::Nbool	# use random schreier method
    extra_options::Ptr{Cvoid}	# void *extra_options;
		# arbitrary extra options
    userautomdata::Ptr{Cvoid}
end

const CONSOLWIDTH = 78
const DISPATCH_GRAPH = Ref{Ptr{Cvoid}}(0)
const ADJACENCIES = Ref{Ptr{Cvoid}}(0)
function __init__()
    DISPATCH_GRAPH[] = cglobal((:dispatch_graph,LIB_FILE),Nothing)
    ADJACENCIES[] = cglobal((:adjacencies,LIB_FILE),Nothing)
    nothing
end
DEFAULTOPTIONS_GRAPH() = optionblk(0,false,false,false,true,false,CONSOLWIDTH,
                    C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,
                    100,0,1,0,DISPATCH_GRAPH[],false,C_NULL,C_NULL)
DEFAULTOPTIONS_DIGRAPH() = optionblk(0,true,false,false,true,false,CONSOLWIDTH,
                    C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,ADJACENCIES[],
                    100,0,999,0,DISPATCH_GRAPH[],false,C_NULL,C_NULL)
DEFAULTOPTIONS(::Type{DenseNautyGraph}) = DEFAULTOPTIONS_GRAPH()
DEFAULTOPTIONS(::Type{DenseNautyDiGraph}) = DEFAULTOPTIONS_DIGRAPH()
DEFAULTOPTIONS(::DenseNautyGraph) = DEFAULTOPTIONS_GRAPH()
DEFAULTOPTIONS(::DenseNautyDiGraph) = DEFAULTOPTIONS_DIGRAPH()

function pprintobject(name, obj, sep=", ")
  print("$name(")
  print(join(map(fn -> "$fn=$(getfield(obj, fn))", fieldnames(typeof(obj))), sep))
  print(")")
end

function Base.show(io::IO, ::MIME"text/plain", options::optionblk)
    pprintobject("optionblk", options)
end

mutable struct statsblk
    grpsize1::Cdouble	# /* size of group is */
    grpsize2::Cint	# /* grpsize1 * 10^grpsize2 */
    numorbits::Cint	# /* number of orbits in group */
    numgenerators::Cint	# /* number of generators found */
    errstatus::Cint	# /* if non-zero : an error code */
    numnodes::Culong	# /* total number of nodes */
    numbadleaves::Culong	# /* number of leaves of no use */
    maxlevel::Cint	# /* maximum depth of search */
    tctotal::Culong	# /* total size of all target cells */
    canupdates::Culong	# /* number of updates of best label */
    invapplics::Culong	# /* number of applications of invarproc */
    invsuccesses::Culong	# /* number of successful uses of invarproc() */
    invarsuclevel::Cint	# /* least level where invarproc worked */
end

const statsblk() = statsblk(zeros(13)...)

function Base.show(io::IO, ::MIME"text/plain", stats::statsblk)
    pprintobject("statsblk", stats)
end

"""
    densenauty(g::DenseNautyGraph, options = DEFAULTOPTIONS_GRAPH(), partition = nothing) -> Tuple\\
    densenauty(g::DenseNautyDiGraph, options = DEFAULTOPTIONS_DIGRAPH(), partition = nothing) -> Tuple

This is the low-level interface to nauty, with a minimal amount of changes to
data. Its arguments match closely the call to densenauty: a graph, an option
block, an optionally a partition, in the form of a pair of lists (lab,ptn):
- `g` is a dense graph or digraph in Nauty format
- `options` is a `struct` following the Nauty structure layout
- `partition` is either `nothing` or `Tuple{Vector{Cint},Vector{Cint}}` representing a partition of the vertices of the graph, in the format `(lab,ptn)`.

The return value is (stats, orbits, ...) a stats block, a description of the
minimal vertices (numbered from 0) in each orbit, and (if options.getcanon) a
list (also numbered from 0) matching the graph's vertices with those in the
canonical graph, and the canonical graph itself.
"""
function densenauty(g::DenseNautyGraph,
               options::optionblk = DEFAULTOPTIONS_GRAPH(),
                    partition::Union{Nothing,Tuple{Vector{Cint},Vector{Cint}}} = nothing)
    __densenauty(g, options, partition)
end
function densenauty(g::DenseNautyDiGraph,
                    options::optionblk = DEFAULTOPTIONS_DIGRAPH(),
                    partition::Union{Nothing,Tuple{Vector{Cint},Vector{Cint}}} = nothing)
    __densenauty(g, options, partition)
end

function __densenauty(g::DenseNautyXGraph, options, partition)
    n = nv(g)
    stats = statsblk()

    if options.defaultptn == 1
        @assert partition == nothing
        lab = zeros(Cint, n)
        ptn = zeros(Cint, n)
    else
        @assert length(partition[1]) == length(partition[2]) == n
        lab = partition[1]
        ptn = partition[2]
    end

    if options.getcanon == 1
        outgraph = typeof(g)(n)
    else
        outgraph = nothing
    end

    orbits = zeros(Cint, n)

    ccall((:densenauty, LIB_FILE), Cvoid,
          (Ptr{UInt64}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ref{optionblk}, Ref{statsblk}, Cint, Cint, Ptr{UInt64}),
          g.data.chunks, lab, ptn, orbits, options, stats, num_setwords(n), n, outgraph == nothing ? C_NULL : outgraph.data.chunks)
    
    if options.getcanon == 1
        result = (lab, outgraph)
    else
        result = ()
    end

    (stats, orbits, result...)
end

list2set(ngroups,l) = IntDisjointSets{Int}(l,zeros(length(l)),ngroups)
                                   
function userautomproc_jl(count::Cint, permptr::Ptr{Cint}, orbitsptr::Ptr{Cint}, numorbits::Cint, stabvertex::Cint, n::Cint, c_data::Ptr{Cvoid})
    data = unsafe_pointer_to_objref(c_data)
    
    perm = unsafe_wrap(Array, permptr, n)
    orbits = unsafe_wrap(Array, orbitsptr, n)

    rec = (permutation = Permutation(perm.+1),orbits = list2set(numorbits,orbits.+1),stabvertex = stabvertex+1)
    push!(data, rec)
    nothing
end

mutable struct NAUTY
    orbits::IntDisjointSets{Int}
    grpsize::BigInt
    generators::Vector{@NamedTuple{permutation::Permutation,orbits::IntDisjointSets{Int},stabvertex::Int}}
    lab::Permutation
    canong::DenseNautyGraph
    NAUTY() = new()
end

function Base.show(io::IO, ::MIME"text/plain", nauty::NAUTY)
    print(io, "(orbits: $(nauty.orbits), group size: $(nauty.grpsize)")
    isdefined(nauty,:lab) && print(io,", labeling: $(nauty.lab)")
    isdefined(nauty,:canong) && print(io,", canong")
    isdefined(nauty,:generators) && print(io,", generators: $([x.permutation for x=nauty.generators])")
    print(io, ")")
end

"""
    nauty(g::DenseNauty[Di]Graph; getcanon = false, automgroup = false, partition = nothing) -> NAUTY

nauty is a higher-level interface to nauty, which relies on densenauty. The
arguments are a graph `g` and optional named arguments getcanon::Bool,
automgroup::Bool and partition, which is either `nothing` or a tuple
`(lab,ptn)` as in the raw Nauty format, or a list of lists of vertices
(numbered from 1).

The return value is a mutable struct NAUTY with entries
- `orbits::IntDisjointSet`, the orbits of the automorphism group numbered from 1
- `grpsize::BigInt`, the group size (possibly with rounding errors)
- if automgroup, `generators::Vector{Tuple{Permutation,IntDisjointSets,Int}}, the generators of the automorphism group, representing each generator with the orbits under the generators up to this one and the stabilized vertex up to now
- if getcanon, `lab::Permutation`, the correspondence between old vertices and new ones, and `canong::typeof(g)`, the canonically labelled graph
"""
function nauty(g::DenseNautyXGraph;
               getcanon::Bool = false,
               automgroup::Bool = false,
               partition::Union{Nothing,Tuple{Vector{Cint},Vector{Cint}},Vector{Vector{Int}}} = nothing)
    options = DEFAULTOPTIONS(g)

    if getcanon
        options.getcanon = true
    end

    if partition == nothing
        labptn = nothing
    elseif isa(partition,Tuple)
        options.defaultptn = false
        labptn = partition
    else
        options.defaultptn = false
        n = nv(g)
        lab = Vector{Cint}(undef,n)
        ptn = ones(Cint,n)
        i = 0
        for p=partition
            @assert !isempty(p)
            lab[i+1:i+length(p)] = p.-1
            i += length(p)
            ptn[i] = 0
        end
        labptn = (lab,ptn)
    end

    result = NAUTY()
    
    if automgroup
        options.userautomproc = @cfunction(userautomproc_jl, Nothing, (Cint, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Ptr{Cvoid}))
        generators = @NamedTuple{permutation::Permutation,orbits::IntDisjointSets{Int},stabvertex::Int}[]
        options.userautomdata = pointer_from_objref(generators)

        rv = densenauty(g, options, labptn)        

        result.generators = generators
    else
        rv = densenauty(g, options, labptn)
    end

    stats = rv[1]
    stats.errstatus == 0 || error("densenauty: error $(stats.errstatus)")
    
    result.orbits = list2set(stats.numorbits,rv[2].+1)
    result.grpsize = div(BigInt(ldexp(significand(stats.grpsize1),precision(Float64)))*BigInt(10)^stats.grpsize2*BigInt(2)^exponent(stats.grpsize1)+BigInt(2)^(precision(Float64)-1),BigInt(2)^precision(Float64))
    if getcanon
        result.lab = Permutation(rv[3].+1)
        result.canong = rv[4]
    end
    result
end

"""NautyIsomorphism() may be supplied as "alg" argument to Graphs.Experimental.has_isomorphism
"""
struct NautyIsomorphism <: Graphs.Experimental.IsomorphismAlgorithm
end

function Graphs.Experimental.has_isomorph(g1::AbstractGraph, g2::AbstractGraph, ::NautyIsomorphism;
                         vertex_relation::Union{Nothing, Function}=nothing,
                         edge_relation::Union{Nothing, Function}=nothing)::Bool
    Graphs.Experimental.could_have_isomorph(g1, g2) && return false

    vertex_relation == nothing || error("didn't code vertex relations yet")
    edge_relation == nothing || error("didn't code edge relations yet")
    
    return nauty(DenseNautyGraph(g1),getcanon=true).canong ==
        nauty(DenseNautyGraph(g2),getcanon=true).canong
end

#add sparse version in the future

end
