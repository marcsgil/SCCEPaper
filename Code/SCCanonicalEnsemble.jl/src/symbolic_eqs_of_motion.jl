function J(d)
    J = zeros(Int,2d,2d)
    for j ∈ axes(J,2), i ∈ axes(J,1)
        if i+d == j
            J[i,j] = -1
        elseif j+d == i
            J[i,j] = 1
        end
    end
    J
end

function get_equations_of_motion(H,d,par_symbol=nothing)
    if isnothing(par_symbol)
        par_symbol = "par"
    end

    Meta.parse( "@variables " * prod("u$n::Real " for n ∈ 1:4d) * "$par_symbol::Real t::Real" ) |> eval

    u = Meta.parse("[" * prod(["u$n, " for n ∈ 1:4d]) * "]") |> eval
    y = Meta.parse("[" * prod(["u$n, " for n ∈ 1:2d]) * "]") |> eval
    x = Meta.parse("[" * prod(["u$n, " for n ∈ 2d+1:4d]) * "]") |> eval
    par = Meta.parse(par_symbol) |> eval

    ℍ = H(x+im*J(d)*y/2,par) + H(x-im*J(d)*y/2,par) |> real

    build_function(vcat(Symbolics.gradient(-ℍ, x),Symbolics.gradient( ℍ, y),), u,par)[2] |> eval
end

function squared_hamiltonian_symbol(H,d)
    Meta.parse( "@variables " * prod("x$n::Real " for n ∈ 1:2d) * "par::Real t::Real" ) |> eval

    x = Meta.parse("[" * prod(["x$n, " for n ∈ 1:2d]) * "]") |> eval
    p = view(x,1:d)
    q = view(x,d+1:2d)

    ham = H(x,par)
    grad_p = Symbolics.gradient(ham, p)
    grad_q = Symbolics.gradient(ham, q)
    lap = sum(n -> Symbolics.derivative(grad_p[n],p[n])*Symbolics.derivative(grad_q[n],q[n]), 1:d)

    build_function(ham^2 - lap/4, x,par,expression=Val{false}) |> eval
end