# Functions for Schubert modules and diagrams
#
# Dave Anderson, July 2024


mutable struct Diagram

  m::Matrix{Int8}

end


function int_to_symb(i::Int8)
    symbols = ['.', "\u25A1",  # "□"
    ]
    return symbols[i+1]
end

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

function Base.size( d::Diagram )

  return size( d.m )

end

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

###
# TO DO write inverse, pos2cols
###

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


function Diagram( w::Vector{Int} )

  return Diagram( Rothe(w) )

end


function random_diagram( a::Int, b::Int )

  m = rand(Bool,a,b)
  m = Int8.(m)

  d = Diagram(m)

  return d
end


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



function hasdescent!( d::Diagram, k::Int )

  if k==size(d.m)[1]
    d.m=vcat(d.m,zeros(Int8,1,size(d.m)[2]))
  end

  b=d.m

  fullinds = filter( j->isfull( findall(x->x==1,b[:,j]), k ), collect(1:size(b)[2]) )

#  fullcols = [ b[:,j] for j in 1:size(b)[2] if isfull( findall(x->x==1,b[:,j]), k ) ]
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

#  nzr = findall( row -> any(row .!= 0), eachrow(d.m) )

#  d.m = d.m[nzr,:]

  d.m = trimd(d.m)

  return true

end


function trimd( m::Matrix{Int8} )

   while iszero( m[end,:] )
     m = m[1:end-1,:]
   end

   nzc = findall( col -> any(col .!= 0), eachcol(m) )
   m = m[:,nzc]

   return m

end


function trimd!( d::Diagram )

  d.m=trimd(d.m)

  return d

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


function row_swap( d::Diagram, k::Int )

  m=d.m

  return Diagram( row_swap(m,k) )


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


function row_swap!( d::Diagram, k::Int )

  m=row_swap!(d.m, k)

  return(d)

end



function skop!( d::Diagram, k::Int )

    if hasdescent!(d,k)
      cols = [ d.m[:,j] for j in 1:size(d)[2] ]
      jj=findfirst( x->(x[k]==1 && x[k+1:end]==zeros(Int8,length(x)-k)), cols )
      d.m[k,jj]=0
    end

    row_swap!(d,k)

end



function skop( d::Diagram, k::Int )

    if hasdescent!(d,k)
      cols = [ d.m[:,j] for j in 1:size(d)[2] ]
      jj=findfirst( x->(x[k]==1 && x[k+1:end]==zeros(Int8,length(x)-k)), cols )
      mm=copy(d.m)
      mm[k,jj]=0

      mm = row_swap(mm,k)
      return Diagram(mm)
    else
      return nothing
    end

end


function rkop( d::Diagram, k::Int )
# without swapping rows

    if hasdescent!(d,k)
      cols = [ d.m[:,j] for j in 1:size(d)[2] ]
      jj=findfirst( x->(x[k]==1 && x[k+1:end]==zeros(Int8,length(x)-k)), cols )
      mm=copy(d.m)
      mm[k,jj]=0

      return Diagram(mm)
    else
      return nothing
    end

end


function isclear( d::Diagram )

  n=size(d)[1]
  all( k -> ( isfull(d,k) || hasdescent!(d,k) ), 1:n )

end


function istransparent( d::Diagram )

  if iszero( d.m ) return true end
  if !isclear( d ) return false end

  for k in descents(d)
    if !istransparent(skop(d,k)) return false end
  end

  return true

end



function mult( d::Diagram, m::Int )
  b=repeat( d.m, inner=(1,m) )

  return Diagram(b)

end



function descents( d::Diagram )

  filter( k->hasdescent!(d,k), collect(1:size(d)[1]) )

end


function character( d::Diagram, R::ZZMPolyRing=xy_ring(size(d)[1]*size(d)[2])[1].ring )

  if iszero( d.m )
    return R(1)
  end

  des = descents( d )

  p1 = R(0)

  x = gens(R)

  for k in des
#    while k<=length(x)
       p1 = p1 + x[k]*Rop( character( skop(d,k), R ), k+1, R )
#    end
  end

  p = p1
  for k in 1:maxvar(p1)
    p1 = Rop( p1, 1, R )
    p = p + p1
  end

  return p
end


function Rop( p::ZZMPolyRingElem, k::Int, R::ZZMPolyRing=parent(p) )

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

  x = gens(R)
  n = length(x)

  if k>n || k<1
    return 0
  end

  p1 = Rop(p,k+1,R) - Rop(p,k,R)

  return( p1/x[k] )
  
end

########
# test some guesses

function check_div_diff( d::Diagram )

#=
  if !isclear(d)
    println("not clear")
    return false
  end
=#

  des = descents(d)

  p=character(d)
  R=parent(p)

  for k in des

    dk=skop(d,k)
    pk=character(dk,R)

    if Rop( pk, k+1, R) != Top( p, k, R )
#    if pk != Rop( ddx(p,k,R), k+1, R )
#    if pk != ddx(p,k,R)
      println(k)
      return false
    end
  end
  return true

end


function check_all_dd( n::Int, m::Int )

  perms = permutations(1:n)
  bv = true

  for w in perms
    d=Diagram(w)
    dm=mult(d,m)
    if !check_div_diff(dm)
      bv=false
      println(w)
    end
  end
  return bv
end


function check_rothe_transparent( n::Int, m::Int )

  perms = permutations(1:n)

  for w in perms
    d=Diagram(w)
    dm = mult(d,m)
    if !istransparent(dm)
      println(w)
      return dm
    end
  end

 return true
end

# not used

#= 
function _fullsort( d::Diagram, k::Int )

  b=d.m

  fullcols = [ b[:,j] for j in 1:size(b)[2] if isfull( findall(x->x==1,b[:,j]), k ) ]
  restcols = setdiff( [ b[:,j] for j in 1:size(b)[2] ], fullcols )

  b=hcat( fullcols..., restcols... )

  return Diagram( b )

end
=#
