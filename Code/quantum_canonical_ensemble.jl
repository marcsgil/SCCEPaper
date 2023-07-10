using SolveSchrodinger

function quantum_energy(θ::Real,Es::AbstractArray)
    sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)
end

function quantum_energy(θs::AbstractArray,Es::AbstractArray)
    map(θ->quantum_energy(θ,Es),θs)
end

function quantum_energy(θ,par,Es::Function)
    quantum_energy(θ,Es(par))
end

function quantum_heat(θ::Real,Es::AbstractArray)
    Z = sum(E->exp(-θ*E),Es)
    U = sum(E->E*exp(-θ*E),Es)/Z
    U2 = sum(E->E^2*exp(-θ*E),Es)/Z
    θ^2*( U2 - U^2 )
end

function quantum_heat(θs::AbstractArray,Es::AbstractArray)
    map(θ->quantum_heat(θ,Es),θs)
end

function quantum_heat(θ,par,Es::Function)
    quantum_heat(θ,Es(par))
end