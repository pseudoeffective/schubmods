# Constructors and methods for Diagrams and Schubert modules
#
# Dave Anderson, July 2024


############################

mutable struct Diagram

  m::Matrix{Int8}

end


######################################
# Extend basic functions for Diagrams
######################################

function int_to_symb(i::Int8)
    symbols = ['.', "\u25A1",  # "□"
    ]
    return symbols[i+1]
end


# extend display for Diagram
function Base.show(io::IO, d::Diagram)
    println(io)
    for i in 1:size(d.m, 1)
        print(" ")
	for j in 1:size(d.m, 2)
            print(io, int_to_symb(d.m[i, j]))
            print(" ")
        end
        println(io)
    end
end


# define size of Diagram
function Base.size( d::Diagram )

  return size( d.m )

end


# extend equality
Base.:(==)( d1::Diagram, d2::Diagram ) = d1.m == d2.m

  


##################
# Constructors
##################

function Diagram( b::Vector{Int8} )
# if there is just one column, given as 01 vector

  m = reshape( b, :, 1 )

  return Diagram( m )

end


function Diagram( pos::Vector{Tuple{Int,Int}} )
# specify the positions of boxes in matrix coordinates

  if isempty(pos)
    return Diagram( zeros(Int8,0,0) )
  end

  a = maximum(first.(pos))
  b = maximum(last.(pos))
  m = zeros(Int8,a,b)

  for (i,j) in pos
    m[i,j]=1
  end

  return Diagram(m)
end


function cols2pos( cols::Vector{Vector{Int}} )
# convert sequence of column sets to positions of boxes

  xys = Tuple{Int,Int}[]
  m = length(cols)

  for j in 1:m
    for i in cols[j]
      push!(xys, (i,j) )
    end
  end

  return xys

end


function Diagram( cols::Vector{Vector{Int}} )
# make diagram from list of columns

  pos = cols2pos( cols )

  return Diagram( pos )

end


function Diagram( b::BPD )
# make diagram from BPD, empty cells to boxes

  xys = findall( x->x==0, b.m )
  pos = [ (xy[1],xy[2]) for xy in xys ]

  return Diagram( pos )

end


function Diagram( dc::Drift )
# make diagram from drift configuration

  xys = findall( x->x==0, dc.m )
  pos = [ (xy[1],xy[2]) for xy in xys ]

  return Diagram( pos )

end


function Diagram( w::Vector{Int} )

  return Diagram( Rothe(w) )

end


function random_diagram( a::Int, b::Int )

  m = rand(Bool,a,b)
  m = Int8.(m)

  d = Diagram(m)

  return d
end


function mult( d::Diagram, r::Int )
# repeat each column of d r times

  b=repeat( d.m, inner=(1,r) )

  return Diagram(b)

end



##########
# Properties
##########

function isfull( v::Vector{Int}, k::Int )

  if (k in v) && !(k+1 in v)
    return false
  else
    return true
  end

end


function isfull( d::Diagram, k::Int )

  b = d.m
  cols = [ b[:,j] for j in 1:size(b)[2] ]

  cols = [ findall( x->x==1, c ) for c in cols ]

  all( x->isfull(x,k), cols )

end  



function is_contained( v1::Vector{Int}, v2::Vector{Int} )

  all( x -> x in v2, v1 )

end


function does_contain( v1::Vector{Int}, v2::Vector{Int} )

  all( x -> x in v1, v2 )

end

function does_strictly_contain( v1::Vector{Int}, v2::Vector{Int} )

  all( x -> x in v1, v2 ) && v1!=v2

end


function hasdescent!( d::Diagram, k::Int )

  if k==size(d.m)[1]
    d.m=vcat(d.m,zeros(Int8,1,size(d.m)[2]))
  end

  b=d.m

  fullinds = filter( j->isfull( findall(x->x==1,b[:,j]), k ), collect(1:size(b)[2]) )

  fullcols = [ b[:,j] for j in fullinds ]
  restcols = [ b[:,j] for j in setdiff( collect(1:size(b)[2]), fullinds) ]

  filter!( !iszero, fullcols )

  restcolsk = [ v[1:k+1] for v in restcols ]
  restcolsk = [ findall( x->x==1, v ) for v in restcolsk ]

  ww = sortperm( restcolsk, lt=does_contain )
  restcolsk = restcolsk[ww]
  restcols=restcols[ww]

  ii=findfirst( x->x[k+1:end]==zeros(Int8,length(x)-k), restcols )

  if ii==nothing
    return false
  else
    v1=restcols[ii]
    v1k=restcolsk[ii]
  end


  restcols = vcat( restcols[1:ii-1], restcols[ii+1:end] )
  restcolsk = vcat( restcolsk[1:ii-1], restcolsk[ii+1:end] )

  for v in restcolsk
    if !does_contain( v1k, v ) || v[end]>k
      return false
    end
  end

  d.m = hcat(fullcols..., v1, restcols... )

  d.m = trimd(d.m)

  return true

end


# rewrite so as not to mutate d
function hasdescent( d::Diagram, k::Int )

  b=d.m

  if k==size(d.m)[1]
    b=vcat(d.m,zeros(Int8,1,size(d.m)[2]))
  end

  fullinds = filter( j->isfull( findall(x->x==1,b[:,j]), k ), collect(1:size(b)[2]) )
  restinds = setdiff( collect(1:size(b)[2]), fullinds)

  # list the full columns
  fullcols = [ b[:,j] for j in fullinds ]
  # list the rest of the columns
  restcols = [ b[:,j] for j in restinds ]

  filter!( !iszero, fullcols )

  # chop off the entries below row k
  restcolsk = [ v[1:k+1] for v in restcols ]
  restcolsk = [ findall( x->x==1, v ) for v in restcolsk ]

  # make sure the entries row k and above are sorted by containment
  ww = sortperm( restcolsk, lt=does_strictly_contain )
  restcolsk = restcolsk[ww]
  restcols = restcols[ww]

  # find a column which stops at row k, for potential border cell
  ii=findfirst( x->x[k+1:end]==zeros(Int8,length(x)-k), restcols )

  if ii==nothing
    return false, (0,0)
  else
    v1=restcols[ii]
    v1k=restcolsk[ii]
  end


  restcols = vcat( restcols[1:ii-1], restcols[ii+1:end] )
  restcolsk = vcat( restcolsk[1:ii-1], restcolsk[ii+1:end] )

  for v in restcolsk
    if !does_contain( v1k, v ) || v[end]>k
      return false, (0,0)
    end
  end

  bc = (k,restinds[ ww[ii] ])
  return true,bc

end



function isclear( d::Diagram )

  n=size(d)[1]
  all( k -> ( isfull(d,k) || hasdescent(d,k)[1] ), 1:n )

end


function istransparent( d::Diagram )

  if iszero( d.m ) return true end
  if !isclear( d ) return false end

  for k in descents(d)
    if !istransparent(skop(d,k)) return false end
  end

  return true

end


function istranslucent( d::Diagram )
# any single-column diagram is translucent

  if iszero( d.m ) return true end
  if size( trimd(d.m) )[2]==1 return true end
  if !isclear( d ) return false end

  for k in descents(d)
    if !istranslucent(skop(d,k)) return false end
  end

  return true

end


function descents( d::Diagram )
# return list of descents of d

  filter( k->hasdescent(d,k)[1], collect(1:size(d)[1]) )

end



function reduced_words( d::Diagram )
# return list of reduced words of d

  if iszero( d.m ) return Vector{Int}[[]] end

  if !isclear(d)
     throw(ArgumentError("d is not transparent"))
  end

  rw = Vector{Int}[]

  for k in descents(d)
    dk = skop(d,k)
    rwk = reduced_words(dk)
    push!.(rwk,k)
    rw = vcat(rw,rwk)
  end

  return rw
end


############
# Operations
############

function trimd( m::Matrix{Int8} )
# remove final empty rows and all empty columns

   while iszero( m[end,:] )
     m = m[1:end-1,:]
   end

   nzc = findall( col -> any(col .!= 0), eachcol(m) )
   m = m[:,nzc]

   return m

end

function trimd( d::Diagram )
# remove final empty rows and all empty columns

  return Diagram( trimd(d.m) )

end


function trimd!( d::Diagram )
# remove final empty rows and all empty columns, in place

  d.m=trimd(d.m)

  return d

end



function row_swap( d::Diagram, k::Int )
# exchange rows k and k+1

  m=d.m

  return Diagram( row_swap(m,k) )


end


function row_swap( m::Matrix{T}, k::Int ) where T

  (a,b)=size(m)

  if k>a || k<1
    return m
  elseif k==a
    mm=vcat(m,zeros(Int8,1,b))
  else
    mm=deepcopy(m)
  end

  rows = vcat( [ mm[i,:] for i in 1:k-1 ], [mm[k+1,:]], [mm[k,:]], [mm[i,:] for i in k+2:a] )

  mm = vcat([r' for r in rows]...)

  return mm

end



function row_swap!( d::Diagram, k::Int )
# exchange rows k and k+1, in place

  row_swap!(d.m, k)

  return(d)

end


function row_swap!( m::Matrix{T}, k::Int ) where T

  (a,b)=size(m)

  if k>a || k<1
    return nothing
  elseif k==a
    m=vcat(m,zeros(T,1,b))
  end

  rows = vcat( [ m[i,:] for i in 1:k-1 ], [m[k+1,:]], [m[k,:]], [m[i,:] for i in k+2:a] )

  m = vcat([r' for r in rows]...)

  return m

end



function skop( d::Diagram, k::Int )
# perform s_k operation on diagram D
    hd,bc=hasdescent(d,k)
    if hd
      (i,j)=bc
      mm=copy(d.m)
      mm[i,j]=0

      mm = row_swap(mm,k)
      return Diagram(mm)
    else
      return nothing
    end

end



function skop!( d::Diagram, k::Int )
# do s_k operation in place

    if hasdescent!(d,k)
      cols = [ d.m[:,j] for j in 1:size(d)[2] ]
      jj=findfirst( x->(x[k]==1 && x[k+1:end]==zeros(Int8,length(x)-k)), cols )
      d.m[k,jj]=0
    end

    row_swap!(d,k)

end



function rkop( d::Diagram, k::Int )
# remove box without swapping rows

    hd,b=hasdescent(d,k)
    if hd
      cols = [ b[:,j] for j in 1:size(b)[2] ]
      jj=findfirst( x->(x[k]==1 && x[k+1:end]==zeros(Int8,length(x)-k)), cols )
      mm=copy(b)
      mm[k,jj]=0

      return Diagram(mm)
    else
      return nothing
    end

end


#############
# Characters of translucent diagrams
#############


function character( d::Diagram, R::ZZMPolyRing=xy_ring(size(d)[1]*size(d)[2])[1].ring )
# character of a translucent diagram, using recursion

  if iszero( d.m )
    return R(1)
  end

  # for single-column diagrams, use single-column method
  if size( trimd(d.m) )[2]==1 
     dmt=reshape( trimd(d.m), : )
     aa = findall( x->x==1, dmt )
     return character(aa,R)
  end
  #

  des = descents( d )

  # check diagram is translucent, otherwise recursion does not work
  notdescent = [ i for i in 1:size(d)[1] if !(i in des) ]
  if !all( k-> isfull(d,k), notdescent )
     throw(ArgumentError("diagram is not translucent"))
  end
  # (note: the above checks d is clear, and will be done at each step of the recursion)

  p1 = R(0)

  x = gens(R)

  # sum over descents
  for k in des
       p1 = p1 + x[k]*Rop( character( skop(d,k), R ), k+1, R )
  end

  # apply 1/(1-R1)
  return Zop(p1,R)
end


function character( aa::Vector{Int}, R::ZZMPolyRing=xy_ring(maximum(aa))[1].ring )
# character of single-column diagram, from basis

  exps = all_bb_less_than(aa)

  p=R(0)
  x=gens(R)
  for bb in exps
    tt=R(1)
    for i in 1:length(bb)
      tt = tt*x[bb[i]]
    end
    p = p + tt
  end

  return p
    
end


function Rop( p::ZZMPolyRingElem, k::Int, R::ZZMPolyRing=parent(p) )
# Bergeron-Sottile operator

  x = gens(R)
  n = length(x)

  if k>n || k<1
    return p
  end

  p1 = evaluate( p, [ x[k] ], [0] )

  for i in k+1:n
    p1 = evaluate( p1, [x[i]], [x[i-1]] )
  end

  return p1

end



function Top( p::ZZMPolyRingElem, k::Int, R::ZZMPolyRing=parent(p) )
# Nadeau-Spink-Tewari trimming operator

  x = gens(R)
  n = length(x)

  if k>n || k<1
    return 0
  end

  p1 = Rop(p,k+1,R) - Rop(p,k,R)

  return( p1/x[k] )
  
end


function Zop( p::ZZMPolyRingElem, R::ZZMPolyRing=parent(p) )
# Nadeau-Spink-Tewari 1/(1-R1) operator

  x = gens(R)
  n = length(x)

  if n<1 return 0 end

  p1 = p

  for k in 1:maxvar(p1)
    p1 = Rop( p1, 1, R )
    p = p + p1
  end

  return p
  
end




function all_bb_less_than(aa::Vector{Int})
    # Check if all elements are positive
    @assert all(a -> a > 0, aa) "All elements must be positive integers"
    # Check if input vector is strictly increasing
    @assert issorted(aa, lt=<) "Input vector must be strictly increasing"
    
    result = Vector{Int}[]
    
    function backtrack(current::Vector{Int}, index::Int, prev::Int)
        if index > length(aa)
            push!(result, copy(current))
            return
        end
        
        for i in (prev + 1):aa[index]
            push!(current, i)
            backtrack(current, index + 1, i)
            pop!(current)
        end
    end
    
    backtrack(Int[], 1, 0)
    
    return result
end

