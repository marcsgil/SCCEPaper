using LinearAlgebra

function normalization(ψs,qs)
    1/√sum(ψ->abs2(ψ)*(qs[2]-qs[1]),ψs) 
end

function get_matrix(Vs,par_p,par_q,Δq)
    N = length(Vs)
    (par_p/Δq^2)*SymTridiagonal(2*ones(N),-ones(N-1)) + par_q*Diagonal(Vs)
end

function solve(V,qs,par_p,par_q)
    Δq = (last(qs)-first(qs))/length(qs)
    Vs = V.(qs)
    Es,ψs = eigen(get_matrix(Vs,par_p,par_q,Δq))
    Ns = [ normalization(view(ψs,:,n),qs) for n in axes(ψs,2) ]
    Es,ψs*Diagonal(Ns)
end