module Banmi

import Mix, Stats

# A struct to store the sorted and indexed data set
#
type DataDescription
    w::Matrix{Int}     # discrete data
    z::Matrix{Float64} # continuous data
    n::Int             # number of rows
    p::Int             # number of discrete variables
    q::Int             # number of continuous variables
    d::Array{Int}      # dimensions of the contingency table
    rw::Matrix{Int}    # missing value indicator for w
    rz::Matrix{Int}    # missing value indicator for z
    nmis::Array{Int}   # number of missing values in each column
end

type Hyperparameters
    dp_weight::Float64

    # Parameters to the unrestricted general location model
    genloc::Mix.UnrestrictedModel

    # Parameters to the priors on the shape variables
    lambda_a::Array{Float64}
    lambda_b::Array{Float64}
    sigma_a::Array{Float64}
    sigma_b::Array{Float64}
end

# A struct to store an instance of the model parameters
#
type Parameters
    x::Matrix{Int}          # discrete modes
    mu::Matrix{Float64}     # continuous means
    lambda::Array{Float64}  # discrete shape parameter
    sigma::Array{Float64}   # continuous shape parameter
end

function ismissing(x::Array{Int})
    x .< 1
end

function ismissing(x::Array{Float64})
    isnan(x)
end

function describe(w::Matrix{Int}, z::Matrix{Float64})
    w = copy(w); z = copy(z)
    n = size(w, 1)
    p = size(w, 2)
    q = size(z, 2)

    # get dimensions of contingency table
    d = mapslices(max, w, 1)

    # missingness indicators for z
    rz = 1 * ismissing(z)
    nmisz = mapslices(sum, rz, 1)
    mdpz = rz * (2 .^ ([1:q] - 1)) + 1
    rz = 1 - rz

    # missingness indicators for w
    rw = 1 * ismissing(w)
    nmisw = mapslices(sum, rw, 1)
    nmis = hcat(nmisw, nmisz)
    mdpw = rw * (2 .^ ([1:p] - 1)) + 1
    rw = 1 - rw

    DataDescription(w, z, n, p, q, d, rw, rz, nmis)
end

function describe(data::Matrix{Float64}, p::Int; mvcode::Float64=0.0)
    w = data[:, 1:p]
    w[w .== mvcode] = 0.0
    w = int(w)

    z = data[:, (p + 1):]
    z[z .== mvcode] = NaN

    describe(w, z)
end

function hyperparameters(desc::DataDescription, 
                         genloc::Mix.UnrestrictedModel, dp_weight::Float64, 
                         lambda_a::Array{Float64}, lambda_b::Array{Float64})

    # set sigma_a and sigma_b to have a mean given by Silverman's rule and a
    # weight equal to the number of complete observations. The 20.0 in the
    # denominator is to downweight the prior.

    sigma_a = Array(Float64, desc.q)
    sigma_b = Array(Float64, desc.q)

    for j in 1:desc.q
        col = desc.z[:, j]
        col = col[!ismissing(col)]
        n_obs = length(col)

        mean_shape = 1.22 * Stats.var(col) * (n_obs ^ -0.4)
        sigma_a[j] = n_obs / 20.0
        sigma_b[j] = sigma_a[j] * mean_shape
    end

    Hyperparameters(dp_weight, genloc, lambda_a, lambda_b, sigma_a, sigma_b)
end

function hyperparameters(desc::DataDescription, 
                         genloc::Mix.UnrestrictedModel,
                         dp_weight::Float64,
                         lambda_a::Float64, lambda_b::Float64)

    la = Array(Float64, desc.p)
    fill!(la, lambda_a)
    lb = Array(Float64, desc.p)
    fill!(lb, lambda_b)

    hyperparameters(desc, genloc, dp_weight, la, lb)
end

# Make an initial draw of parameters
#
function draw_parameters(desc::DataDescription, hp::Hyperparameters)
    x = Array(Int, desc.n, desc.p)
    mu = Array(Float64, desc.n, desc.q)

    for i in 1:desc.n
        # draw mean from hp.genloc
        # draw multivariate normal
    end
end

# Initialize the missing values by drawing them from the Mix model
#
function init_missing_values(w::Matrix{Int}, z::Matrix{Float64},
                             desc::DataDescription, theta::Parameters)

    wout = copy(w); zout = copy(z)                         

    # first fill in the categorical values

    wout, zout
end

# Update the values in theta
#
function data_augmentation!(theta::Parameters, desc::DataDescription, steps::Int)
end

function impute(data, desc, theta)
end

end # module
