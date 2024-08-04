# schubmods
*julia code for Schubert and diagram modules

This file provides methods for computing characters of Schur modules for diagrams.  The polynomials arising as such characters include the [Schubert polynomials](https://en.wikipedia.org/wiki/Schubert_polynomial) in the case the given diagram is the [Rothe diagram](https://en.wikipedia.org/wiki/Permutation#Numbering_permutations) of a permutation, and in particular [Schur polynomials](https://en.wikipedia.org/wiki/Schur_polynomial) as special cases.  

To use it, make sure the file is where Julia can find it, say in your current working directory, and load the [SchubertPolynomials](https://github.com/pseudoeffective/SchubertPolynomials.jl/blob/main/README.md) package:

```julia-repl
julia> using SchubertPolynomials
julia> include("schubmods.jl")
```
Examples of diagrams and characters will be given below.  (Soon!)

