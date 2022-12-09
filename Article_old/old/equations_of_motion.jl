using Symbolics

function equations_of_motion(hamiltonian)
    @syms p::Real q::Real yp::Real yq::Real par::Real
    args = [yp,yq,p,q,par]
    Dyp = Differential(yp)
    Dyq = Differential(yq)
    Dp = Differential(p)
    Dq = Differential(q)

    ℌ = hamiltonian([p-im*yq/2,q+im*yp/2],par) + hamiltonian([p+im*yq/2,q-im*yp/2],par)

    dyp = build_function( Dyp(ℌ) |> expand_derivatives |> real, args) |> eval
    dyq = build_function( Dyq(ℌ) |> expand_derivatives |> real, args) |> eval
    dp = build_function( Dp(ℌ) |> expand_derivatives |> real, args) |> eval
    dq = build_function( Dq(ℌ) |> expand_derivatives |> real, args) |> eval

    fy(y,x,par) = [-dp([y[1],y[2],x[1],x[2],par]),-dq([y[1],y[2],x[1],x[2],par])]
    fx(y,x,par) = [dyp([y[1],y[2],x[1],x[2],par]),dyq([y[1],y[2],x[1],x[2],par])]
    fy,fx
end