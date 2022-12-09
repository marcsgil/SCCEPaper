function symetric_integer_range(N)
    #=Returns a range of integers of length N that is symetric with respect to 0 when N is odd, and has an extra negative number when N is even.

    Examples:
    symetric_integer_range(4) = -2:1
    symetric_integer_range(5) = -2:2
    =#

    -N÷2:N - N÷2 - 1
end

struct DFTGrid{ndims<:Val}
    N::Int64
    max_val::Float64
    DFTGrid(N,max_val,ndims) = new{Val{ndims}}(N,max_val)
end

interval(grid::DFTGrid) = grid.max_val/(grid.N÷2)
reciprocal_interval(grid::DFTGrid) = π/grid.max_val
direct_grid(grid::DFTGrid) = interval(grid)*symetric_integer_range(grid.N)
reciprocal_grid(grid::DFTGrid) = reciprocal_interval(grid)*symetric_integer_range(grid.N)

function multi_dim_direct_grid(grid::DFTGrid{Val{ndims}}) where ndims
    freqs = fill(ntuple(j->zero(typeof(grid.max_val)), ndims), ntuple(j->grid.N,ndims))

    cart_indices = CartesianIndices(freqs)
    dir_grid = direct_grid(grid)

    for n in eachindex(freqs)
        freqs[n] = ntuple(j->dir_grid[cart_indices[n][j]], ndims) 
    end

    freqs
end

function multi_dim_reciprocal_grid(grid::DFTGrid{Val{ndims}}) where ndims
    freqs = fill(ntuple(j->zero(typeof(grid.max_val)), ndims), ntuple(j->grid.N,ndims))

    cart_indices = CartesianIndices(freqs)
    rec_grid = reciprocal_grid(grid)

    for n in eachindex(freqs)
        freqs[n] = ntuple(j->rec_grid[cart_indices[n][j]], ndims) 
    end

    freqs
end

function evaluate_at_direct_grid(f::Function,grid::DFTGrid;perform_shift=false,use_GPU=false)
    if use_GPU
        if perform_shift
            ifftshift(map(x->f(x),CuArray(multi_dim_direct_grid(grid))))
        else
            map(x->f(x),CuArray(multi_dim_direct_grid(grid)))
        end
    else
        if perform_shift
            ifftshift(map(x->f(x),multi_dim_direct_grid(grid)))
        else
            map(x->f(x),multi_dim_direct_grid(grid))
        end
    end
end

function evaluate_at_reciprocal_grid(f::Function,grid::DFTGrid;perform_shift=false,use_GPU=false)
    if use_GPU
        if perform_shift
            ifftshift(map(x->f(x),CuArray(multi_dim_reciprocal_grid(grid))))
        else
            map(x->f(x),CuArray(multi_dim_reciprocal_grid(grid)))
        end
    else
        if perform_shift
            ifftshift(map(x->f(x),multi_dim_reciprocal_grid(grid)))
        else
            map(x->f(x),multi_dim_reciprocal_grid(grid))
        end
    end
end