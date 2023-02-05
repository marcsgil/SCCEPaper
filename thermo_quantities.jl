function energy_output(sol,i,θ,par,X,H)
    Z = exp( sol[end].Δ - θ*H(view(X,:,i),par) )*√abs(det(sol[end].jac_x))
    [Z,Z*H(sol[end].x,par)],false
end

function energy(θ::Number,par,fy,fx,H)
    prob = IntegralProblem((X,par)->integrand(X,θ,par,fy,fx,(sol,i)->energy_output(sol,i,θ,par,X,H)),fill(-Inf,2),fill(Inf,2),par,nout=2)
    sol = solve(prob,CubatureJLh(),reltol=1e-3,abstol=1e-4)
    sol[2]/sol[1]
end

energy(θ::Number,par,H) = energy(θ::Number,par,equations_of_motion(H)...,H)

function energy(θs::AbstractArray,par,fy,fx,H)
    ThreadsX.map(θ->energy(θ::Number,par,fy,fx,H),θs)
end

energy(θ::AbstractArray,par,H) = energy(θ::AbstractArray,par,equations_of_motion(H)...,H)

function heat_output(sol,i,θ,par,X,H,H2)
    Z = exp( sol[end].Δ - θ*H(view(X,:,i),par) )*√abs(det(sol[end].jac_x))
    [Z,Z*H(sol[end].x,par),Z*H2(sol[end].x,par)],false
end

function heat(θ::Number,par,fy,fx,H,H2)
    prob = IntegralProblem((X,par)->integrand(X,θ,par,fy,fx,(sol,i)->heat_output(sol,i,θ,par,X,H,H2)),-20*ones(2),20*ones(2),par,batch=16,nout=3)
    sol = solve(prob,CubatureJLh(),reltol=1e-5,abstol=1e-5)
    θ^2*(sol[3]/sol[1]-(sol[2]/sol[1])^2)
end

function heat(θs::AbstractArray,par,fy,fx,H,H2)
    ThreadsX.map(θ->heat(θ::Number,par,fy,fx,H,H2),θs)
end