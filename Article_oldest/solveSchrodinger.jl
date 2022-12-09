using LinearAlgebra,SparseArrays,Arpack,Kronecker

D(N) = SymTridiagonal(-2ones(N),ones(N-1))

function solveSchrodinger(xs,V;mass=1,par=nothing)
    H = -D(length(xs))/(2*mass*(xs[2]-xs[1])^2) + diagm( map(x->V(x,par), xs) )
    eigen(H)
end

function solveSchrodinger(xs,ys,V;xmass=1,ymass=1,nev=50,par=nothing)
    H = sparse( -( D(length(xs))/(2*xmass*(xs[2]-xs[1])^2) ⊕ D(length(ys))/(2*ymass*(ys[2]-ys[1])^2) )
    + diagm(vec(map(r->V(r...,par), Iterators.product(xs,ys)))) )
    eigs(H,nev=nev,which=:SR)
end

#=
Harmonic oscilator test

V(x,ω) = ω*x^2/2
V(x,y,ω) = ω*(x^2+y^2)/2

N=128
xs = LinRange(-5,5,N)
ys = LinRange(-5,5,N)
Es,ψs = solveSchrodinger(xs,ys,V,par=1)
=#