using DiffEqGPU, DifferentialEquations, StaticArrays, Symbolics, Latexify

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

function get_equations_of_motion(H,d,index)
    Meta.parse( "@variables " * prod("u$n::Real " for n in 1:4d) * "par::Real" ) |> eval

    exp_y = Meta.parse("[" * prod(["u$n, " for n in 1:2d]) * "]")
    exp_x = Meta.parse("[" * prod(["u$n, " for n in 2d+1:4d]) * "]")

    y = eval(exp_y)
    x = eval(exp_x)

    ℍ = H(x+im*J(d)*y/2,par) + H(x-im*J(d)*y/2,par) |> real

    if index ≤ 2d
        Symbolics.derivative(-ℍ, eval(Meta.parse("u$(2d + index)"))) |> string
    else
        Symbolics.derivative(ℍ,  eval(Meta.parse("u$(index-2d)"))) |> string
    end
end
##
H(x,ω) = ω*sum(x->x^2,x)/2

get_equations_of_motion(H,1,4)

function translate(expr)
    replace(expr,"u1"=>"u[1]","u2"=>"u[2]","u3"=>"u[3]","u4"=>"u[4]","par"=>"par[1]")
end


Meta.parse(get_equations_of_motion(H,1,4) |> translate) |> eval
##
function f(du, u, par, t)
    Meta.parse("du[1] = " * get_equations_of_motion(H,1,1) |> translate) |> eval
    #du[2] = Meta.parse(get_equations_of_motion(H,1,2) |> translate) |> eval
    #du[3] = Meta.parse(get_equations_of_motion(H,1,3) |> translate) |> eval
    #du[4] = Meta.parse(get_equations_of_motion(H,1,4) |> translate) |> eval
end

f([0,0,0,0], [0,0,0,0], [0], 0)

u0 = rand(Float32,4)
tspan = (0f0,1f0)
ω = [1f0]
prob = ODEProblem(f,u0,tspan,ω)
prob_func = (prob,i,repeat) -> remake(prob, u0 = rand(Float32,4))
monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy=false)
sol = solve(monteprob,Tsit5(),trajectories=10_000,saveat=1.0f0)

