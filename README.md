# schubmods
*julia code for Schubert and diagram modules*

This file provides methods for computing characters of Schur modules for diagrams.  The polynomials arising as such characters include the [Schubert polynomials](https://en.wikipedia.org/wiki/Schubert_polynomial) in the case the given diagram is the [Rothe diagram](https://en.wikipedia.org/wiki/Permutation#Numbering_permutations) of a permutation, and in particular [Schur polynomials](https://en.wikipedia.org/wiki/Schur_polynomial) as special cases.  

To use it, make sure the file is where Julia can find it, say in your current working directory, and load the [SchubertPolynomials](https://github.com/pseudoeffective/SchubertPolynomials.jl/blob/main/README.md) package:

```julia-repl
julia> using SchubertPolynomials
julia> include("schubmods.jl")
```

## Examples

There are several ways to contruct a diagram.  A basic one is by specifying positions (in matrix coordinates).
```julia-repl
julia> pos = [ (1,1), (1,2), (2,3) ];

julia> d = Diagram(pos)

 □ □ . 
 . . □ 
```
One can also get the Rothe diagram of a permutation.
```julia-repl
julia> w = [1,4,6,2,5,3];

julia> dw = Diagram(w)

 . . . . . 
 . □ □ . . 
 . □ □ . □ 
 . . . . . 
 . . □ . . 
```
The function `mult(d,m)` repeats each column of `d`.
```julia-repl
julia> mult(dw,3)

 . . . . . . . . . . . . . . . 
 . . . □ □ □ □ □ □ . . . . . . 
 . . . □ □ □ □ □ □ . . . □ □ □ 
 . . . . . . . . . . . . . . . 
 . . . . . . □ □ □ . . . . . . 
```

The *descents* of a diagram `d` are the row indices where `d` satisfies a certain property.  The definition is motivated by requiring `Diagram(w)` to have descents equal to the (right) descents of `w`.
```julia-repl
julia> descents(dw)
2-element Vector{Int64}:
 3
 5
```
*Warning: the internal function `hasdescents!` mutates the diagram, so `d` may look different after calling `descents`.  The intended applications are insensitive to permutations of columns, though.*

If `k` is a descent of `d`, the function `skop(d,k)` removes a certain box in row `k` and swaps rows `k` and `k+1`.  When `d` is the Rothe diagram of `w`, the result is the Rothe diagram of `ws_k`.
```julia-repl
julia> dw3=skop(dw,3)

 . . . 
 □ □ . 
 . . . 
 . □ □ 
 . □ . 


julia> w3=[1,4,2,6,5,3];

julia> Diagram(w3)

 . . . . . 
 . □ □ . . 
 . . . . . 
 . . □ . □ 
 . . □ . . 
```
The function `trimd` removes empty columns.
```julia-repl
julia> trimd(Diagram(w3))==dw3
true
```
It can happen that `k` is still a descent of `skop(d,k)` (though not for Rothe diagrams).
```
julia> d

 □ □ . 
 . . □ 

julia> descents(d)
2-element Vector{Int64}:
 1
 2

julia> d1=skop(d,1)

 □ . . 
 . . □ 

julia> d11=skop(d1,1)

 □ . 
 . . 
```

A diagram is *translucent* if it satisfies a certain recursively defined property; translucent diagrams include all Rothe diagrams of permutations, as well as some more.  For translucent diagrams, the function `character` computes the character of the corresponding flagged Schur module.
```julia-repl
julia> istranslucent(d)
true

julia> character(d)
x1^3 + x1^2*x2
```
If `d=Diagram(w)` is a Rothe diagram, a theorem of Kraskiewicz and Pragacz says its character is equal to a Schubert polynomial.  The ambient polynomial ring can be specified, which is useful when comparing polynomials obtained by different methods.
```julia-repl
julia> R,x,y=xy_ring(6);
julia> Rx=R.ring;
julia> pw = character(dw,Rx);
julia> spw=schub_poly(w,R);
julia> pw==spw
true
```
