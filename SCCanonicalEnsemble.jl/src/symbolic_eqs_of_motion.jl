function J(d)
    J = zeros(Int,2d,2d)
    for j in axes(J,2)
        for i in axes(J,1)
            if i+d == j
                J[i,j] = -1
            elseif j+d == i
                J[i,j] = 1
            end
        end
    end
    J
end

function get_equations_of_motion(H,d,par_symbol=nothing)
    
    if isnothing(par_symbol)
        par_symbol = "par"
    end

    Meta.parse( "@variables " * prod("u$n::Real " for n in 1:4d) * "$par_symbol::Real t::Real" ) |> eval

    u = Meta.parse("[" * prod(["u$n, " for n in 1:4d]) * "]") |> eval
    y = Meta.parse("[" * prod(["u$n, " for n in 1:2d]) * "]") |> eval
    x = Meta.parse("[" * prod(["u$n, " for n in 2d+1:4d]) * "]") |> eval
    par = Meta.parse(par_symbol) |> eval

    ℍ = H(x+im*J(d)*y/2,par) + H(x-im*J(d)*y/2,par) |> real

    build_function(vcat(Symbolics.gradient(-ℍ, x),Symbolics.gradient( ℍ, y),), u,par)[2] |> eval
end