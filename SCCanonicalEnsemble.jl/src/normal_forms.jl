#Symbolic part
function get_converted_monomials(N)
    # monomials[n] =  ∑ cⱼ (p²+q²)ʲ, the Wigner representation of (̂p²+̂q²)ⁿ
    p = Polynomial([0,1])
    monomials = [p for n in 1:N+1]
    monomials[1] = one(p)
    for n in 3:N+1
        monomials[n] = p*monomials[n-1]-derivative(monomials[n-1])-p*derivative(derivative(monomials[n-1]))
    end
    monomials
end

function normalize_and_change_basis(p::Polynomial)
    #Given a basis 1,…,xⁿ and a polynomial p(x) = ∑ cₙ xⁿ, returns returns p/(2^degrre(p)) with respect to the basis 1,…,uⁿ, where u = x/2

    #This is useful because it is straightforward to obtain the Wigner representation of (̂p²+̂q²)ⁿ, which is done by get_converted_monomials.
    #On the other hand, to perform other calculations, it is easier to express normal forms as powers of [(̂p²+̂q²)/2]ⁿ or [(p²+q²)/2]ⁿ

    ds = coeffs(p)/1
    n = degree(p)
    for j ∈ eachindex(ds)
        ds[j]/=2^(n-j+1)
    end
    Polynomial(ds)
end

function wigner_symbol_of_normal_form(quantum_nf::Polynomial;assume_divided_by_two_basis=true)
    #Returns the Wigner symbol of a quantum normal form
    #If assume_divided_by_two_basis=true, we have quantum_nf = ∑ cₙ ̂uⁿ and result = ∑ dₙ uⁿ, where ̂u = (̂p²+̂q²)/2 and u = (p²+q²)/2
    #If assume_divided_by_two_basis=false, we have quantum_nf = ∑ cₙ ̂uⁿ and result = ∑ dₙ uⁿ, where ̂u = p²+̂q² and u = p²+q²

    converted_monomials = get_converted_monomials(degree(quantum_nf))

    if assume_divided_by_two_basis
        result = zero(quantum_nf)/1
        for (n,c) in enumerate(coeffs(quantum_nf))
            result += c*normalize_and_change_basis(converted_monomials[n])
        end
    else
        result = zero(quantum_nf)
        for (n,c) in enumerate(coeffs(quantum_nf))
            result += c*converted_monomials[n]
        end
    end

    result
end

#Numeric part

function Z_integrand_NF(J,θ,F::Polynomial,ω::Polynomial,ω′::Polynomial)
    #Calculates the integrand of the partition function
    ϕ = θ*ω(J)
    half_ϕ = ϕ/2
    Sᴱ = (ϕ-sinh(ϕ))*J - θ*F(J)
    modified_cosh = ( 1+exp(-ϕ) )/2
    sqrt_term = √abs((1+J*ω′(J)*θ*tanh(half_ϕ)))
    sqrt_term*modified_cosh*exp(Sᴱ + half_ϕ)
end

#=function Z_NF(θ,F::Polynomial,ω::Polynomial,ω′::Polynomial)
    prob = IntegralProblem((J,par)->Z_integrand_NF(J,θ,F::Polynomial,ω::Polynomial,ω′::Polynomial),
            0,100,)
    solve(prob,QuadGKJL()).u
end

function energy_NF(θ,F::Polynomial,ω::Polynomial,ω′::Polynomial)
    -ForwardDiff.derivative(θ->log(Z_NF(θ,F,ω,ω′)),θ)
end=#

function A_integrand(J,θ,F::Polynomial,ω::Polynomial,ω′::Polynomial,G::Polynomial)
    #Calculates the integrand of the expectation value of an operator A = G(J)
    ϕ = θ*ω(J)
    half_ϕ = ϕ/2
    Sᴱ = (ϕ-sinh(ϕ))*J - θ*F(J)
    modified_cosh = ( 1+exp(-ϕ) )/2
    sqrt_term = √abs((1+J*ω′(J)*θ*tanh(half_ϕ)))
    sqrt_term*sum( pair-> pair[2]*J^(pair[1])*modified_cosh^(2*pair[1]+1)*exp(Sᴱ + (2*pair[1]+1)*half_ϕ), pairs(G))
end

function ode_step(y,(θ,F,ω,ω′,G),J)
    dy1 = Z_integrand_NF(J,θ,F,ω,ω′)
    dy2 = A_integrand(J,θ,F,ω,ω′,G)
    SA[dy1,dy2]
end

function expectation_value_NF( θ,F::Polynomial,G::Polynomial,J_max=200.)
    ω = derivative(F)
    ω′ = derivative(ω)
    prob = ODEProblem(ode_step,SA[zero(J_max),zero(J_max)],(0,J_max),(θ,F,ω,ω′,G))
    sol = solve(prob,save_everystep=false,save_start=false)
    sol[2,end]/sol[1,end]
end

function expectation_value_NF( θs::AbstractArray,F::Polynomial,G::Polynomial,J_max=200.)
    ω = derivative(F)
    ω′ = derivative(ω)
    prob = ODEProblem(ode_step,zeros(SVector{2,typeof(J_max)}),(0,J_max),(θs[1],F,ω,ω′,G))

    function prob_func(prob,i,repeat)
        remake(prob,p=(θs[i],F,ω,ω′,G))
    end

    output_func(sol,i) = (sol[2,end]/sol[1,end],false)

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)
    solve(ensemble_prob,trajectories=length(θs))
end

function energy_NF(θs::AbstractArray,H_operator::Polynomial,J_max=200.)
    F = wigner_symbol_of_normal_form(H_operator)
    expectation_value_NF( θs,F,F,J_max)
end

function heat_capacity_NF(θs::AbstractArray,H_operator::Polynomial,J_max=200.)
    F = wigner_symbol_of_normal_form(H_operator)
    F2 = wigner_symbol_of_normal_form(H_operator^2)
    Us = expectation_value_NF( θs,F,F,J_max)
    Us2 = expectation_value_NF( θs,F,F2,J_max)
    map( (U2,U,θ)-> θ^2*(U2-U^2), Us2,Us,θs )
end