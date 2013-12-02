module Mix

import Base.bsxfun

function rngseed(seed::Int)
    if seed <= 0
        throw(DomainError())
    end

    ccall((:rngs_, "libmix"), Void, (Ptr{Cint},), &int32(seed))
end

function na_to_code!{T <: Integer}(x::Array{T}, code::T)
    x[x .< 1] = code
    x
end

function na_to_code!{T <: Real}(x::Array{T}, code::T)
    x[isnan(x)] = code
    x
end

function duplicated{T}(x::Array{T, 1})
    result = falses(length(x))
    seen = Set{T}()
    for (i, v) in enumerate(x)
        if in(v, seen)
            result[i] = true
        else
            push!(seen, v)
        end
    end

    result
end

type MixDataDescription
    w::Matrix{Int32}
    n::Int32
    p::Int32
    d::Array{Int32}
    jmp::Array{Int32}
    z::Matrix{Float64}
    q::Int32
    r::Matrix{Int32}
    rz::Matrix{Int32}
    rw::Matrix{Int32}
    nmis::Array{Int32}
    ro::Array{Int32}
    mdpzgrp::Array{Int32}
    mdpwgrp::Array{Int32}
    mobs::Array{Int32}
    mobsst::Array{Int32}
    nmobs::Array{Int32}
    ncells::Int32
    ngrp::Int32
    npattz::Int32
    npattw::Int32
    npsi::Int32
    psi::Matrix{Int32}
    xbar::Array{Float64}
    sdv::Array{Float64}
end

function describe(w::Matrix{Int}, z::Matrix{Float64})
    w = copy(w); z = copy(z)
    n = size(w, 1)
    p = size(w, 2)
    q = size(z, 2)

    # get dimensions of contingency table
    d = mapslices(maximum, w, 1)

    # missingness indicators for z
    rz = 1 * isnan(z)
    nmisz = mapslices(sum, rz, 1)
    mdpz = rz * (2.^([1:q] - 1)) + 1
    rz = 1 - rz

    # get missing data patterns for w
    rw = 1 * (w .< 1)
    nmisw = mapslices(sum, rw, 1)
    nmis = hcat(nmisw, nmisz)
    mdpw = rw * (2.^([1:p] - 1)) + 1
    rw = 1 - rw

    # calculate the known part of the cell number for rows of w
    cumd = cumprod(d, 2)
    mobs = 1 + ((w .- 1) .* rw) * div(cumd, d)'

    # do row sort
    ro = sortrows(hcat(mdpz, mdpw, mobs, [1:n]))[:, 4]
    w = w[ro, :]
    mdpw = mdpw[ro]
    mobs = mobs[ro]
    rw = rw[ro, :]
    z = z[ro, :]
    mdpz = mdpz[ro]
    rz = rz[ro, :]
    ro = sortperm(ro)

    # compress missing data patterns
    mdpzst = (1:n)[!duplicated(mdpz)]
    mdpz = mdpz[!duplicated(mdpz)]
    npattz = length(mdpz)
    mdpzfin =  push!(mdpzst[2:npattz] - 1, n)
    mdpzgrp = Int[]; mdpwst = Int[]; mdpwtmp = Int[]

    for i in 1:npattz
        tmp = mdpw[mdpzst[i]:mdpzfin[i]]
        mdpwst = vcat(mdpwst,
                      (1:length(tmp))[!duplicated(tmp)] + mdpzst[i] - 1)
        tmp = tmp[!duplicated(tmp)]
        push!(mdpzgrp, length(tmp))
        mdpwtmp = vcat(mdpwtmp, tmp)
    end

    mdpw = mdpwtmp
    npattw = length(mdpw)
    mdpwfin = push!(mdpwst[2:npattw] - 1, n)
    mdpwgrp = Int[]; mobsst = Int[]; mobstmp = Int[]

    for i in 1:npattw
        tmp = mobs[mdpwst[i]:mdpwfin[i]]
        mobsst = vcat(mobsst,
                      (1:length(tmp))[!duplicated(tmp)] + mdpwst[i] - 1)
        tmp = tmp[!duplicated(tmp)]
        push!(mdpwgrp, length(tmp))
        mobstmp = vcat(mobstmp, tmp)
    end

    mobs = mobstmp
    ngrp = length(mobs)

    # create r-matrix for display purposes
    r = hcat(rw, rz)[mdpwst, :]
    ncells = cumd[p]
    jmp = int(trunc(cumd ./ d))
    rz = rz[mdpzst, :]
    rw = rw[mdpwst, :]

    # form matrix of packed storage indices
    npsi = div(q * (q + 1), 2)
    psi = Array(Cint, q, q)
    ccall((:mkpsi_, "libmix"), Void, (Ptr{Cint}, Ptr{Cint}),
          &convert(Cint, q-1), psi)

    # center and scale the columns of z
    mvcode = maximum(z[!isnan(z)]) + 1000
    na_to_code!(z, mvcode)
    xbar = Array(Cdouble, q)
    sdv = Array(Cdouble, q)
    ccall((:ctrsc_, "libmix"), Void,
          (Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          z, &n, &q, xbar, sdv, &mvcode)
    z[z .== mvcode] = NaN

    # return MixDataDescription struct
    nmobs = vcat(mobsst[2:], n + 1) - mobsst
    MixDataDescription(
        int32(w), int32(n), int32(p), int32(d), int32(jmp), z, int32(q),
        int32(r), int32(rz), int32(rw), int32(nmis), int32(ro),
        int32(mdpzgrp), int32(mdpwgrp),
        int32(mobs), int32(mobsst), int32(nmobs), int32(ncells), int32(ngrp),
        int32(npattz), int32(npattw), int32(npsi), int32(psi),
        xbar, sdv)
end

function describe(data::Matrix{Float64}, p::Int; mvcode=0.0)
    w = data[:, 1:p]
    w[w .== mvcode] = 0.0
    w = int(w)

    z = data[:, (p + 1):]
    z[z .== mvcode] = NaN

    describe(w, z)
end

abstract MixModel

type UnrestrictedModel <: MixModel
    sigma::Array{Float64}
    mu::Array{Float64}
    pi::Array{Float64}
end

# TODO implement version that takes a 'start' parameter
#
function em(desc::MixDataDescription; 
            prior=1, maxits::Int=1000, showits::Bool=true, epsilon::Float64=0.0001)

    if length(prior) == 1
        tmp = prior
        prior = Array(Float64, desc.ncells)
        fill!(prior, tmp)
    end

    prior = float64(prior)
    w = !isnan(prior)
    prior[!w] = -999.0
    z = copy(desc.z); na_to_code!(z, 999.0)
    tp = zeros(Int32, desc.p)
    tq = zeros(Int32, desc.q)
    kn1 = zeros(Float64, desc.npsi)
    kn2 = zeros(Float64, desc.q, desc.ncells)
    kn3 = zeros(Float64, desc.ncells)

    ccall((:tobsm_, "libmix"), Void,
         # q          psi        npsi       t1            ncells     t2            t3            npattz
          (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
         # rz         mdpzgrp    npattw     p          rw         mdpwgrp    ngrp       mobs       mobsst
           Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
         # nmobs      n          z             ocw        ocz
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
           &desc.q, desc.psi, &desc.npsi, kn1, &desc.ncells, kn2, kn3,
           &desc.npattz, desc.rz, desc.mdpzgrp, &desc.npattw, &desc.p, desc.rw, desc.mdpwgrp,
           &desc.ngrp, desc.mobs, desc.mobsst, desc.nmobs, &desc.n, z, copy(tp), copy(tq))

    # this block replaces a check for the 'start' parameter
    sigma = zeros(Float64, desc.npsi)
    mu = zeros(Float64, desc.q, desc.ncells)
    ccall((:stvlm_, "libmix"), Void,
          (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
          &desc.q, desc.psi, &desc.npsi, sigma, &desc.ncells, mu)
    pii = ones(Float64, desc.ncells)
    pii[!w] = 0.0

    converged = false
    it = 0
    t1 = copy(sigma)
    t2 = copy(mu)
    t3 = copy(pii)
    if showits; println("Steps of EM:"); end

    while !converged && (it < maxits)
        it += 1
        if showits; print(it, "..."); end
        if it > 1
            sigma[:] = t1
            mu[:] = t2
            pii[:] = t3
        end

        # The following variables are overwritten by the call to estepm, send
        # copies instead: tp, tq, second appearance of kn3, pii

        ccall((:estepm_, "libmix"), Void,
             # q          psi        npsi       ncells     sigma         mu            pii           kn1           kn2           kn3
              (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
             # t1            t2            t3            npattz     rz         mcz        ocz        mdpzgrp    npattw     p          rw         mcw
               Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
             # mdpwgrp    ngrp       mobs       mobsst     nmobs      n          z             d          jmp        c          theta
               Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
               &desc.q, desc.psi, &desc.npsi, &desc.ncells, sigma, mu, copy(pii), kn1,
               kn2, kn3, t1, t2, t3, &desc.npattz, desc.rz, copy(tq), copy(tq), desc.mdpzgrp, &desc.npattw,
               &desc.p, desc.rw, copy(tp), desc.mdpwgrp, &desc.ngrp,
               desc.mobs, desc.mobsst, desc.nmobs, &desc.n, z, desc.d, desc.jmp, copy(tp), copy(kn3))

        ccall((:mstepm_, "libmix"), Void,
              (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
              &desc.q, desc.psi, &desc.npsi, &desc.ncells, t1, t2, t3, &desc.n, prior)

        if any(t3 .< 0)
            error("Estimate outside the parameter space. Check prior.")
        end

        t3[w & (t3 .< 0.0000001/sum(w))] = 0
        c3 = all(abs(t3 - pii) .<= epsilon * abs(pii))
        c2 = all(abs(t2 - mu) .<= epsilon * abs(mu))
        c1 = all(abs(t1 - sigma) .<= epsilon * abs(sigma))
        converged = c1 && c2 && c3
    end

    if showits; println(); end
    UnrestrictedModel(t1, t2, t3)
end

# TODO optinally return correlation matrix?
# 
function getparam(s, theta)
    pii = reshape(theta.pi, int(s.d)...)
    mu = broadcast(.+, broadcast(.*, theta.mu, s.sdv), s.xbar)
    sigma = theta.sigma[s.psi]
    tmp = repmat(s.sdv, 1, int(s.q))
    sigma = sigma .* tmp .* tmp'
    mu[:, reshape(pii .== 0, prod(s.d))] = NaN

    UnrestrictedModel(sigma, mu, pii)
end

function impute(desc::MixDataDescription, theta::MixModel, data; mvcode=0.0)
    sigma = copy(theta.sigma)
    mu = copy(theta.mu)
    pii = copy(theta.pi)
    z = copy(desc.z); na_to_code!(z, 999.0)
    w = copy(desc.w); na_to_code!(w, int32(999))
    tp = zeros(desc.p)
    tq = zeros(desc.q)

    ccall((:istepm_, "libmix"), Void,
         # q          psi        npsi       ncells     sigma         mu            pii           kn1           kn2           kn3
          (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
         # t1            t2            t3            npattz     rz         mcz        ocz        mdpzgrp    npattw     p          rw         mcw 
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
         # mdpwgrp    ngrp       mobs       mobsst     nmobs      n          z             d          jmp        c          theta         chf           w          zz
           Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
           &desc.q, desc.psi, &desc.npsi, &desc.ncells, sigma, mu, pii, 
           copy(sigma), copy(mu), copy(pii), copy(sigma), copy(mu), copy(pii), &desc.npattz, desc.rz, tq, copy(tq), desc.mdpzgrp, 
           &desc.npattw, &desc.p, desc.rw, tp, desc.mdpwgrp, &desc.ngrp,
           desc.mobs, desc.mobsst, desc.nmobs, &desc.n, z, desc.d, desc.jmp, copy(tp), copy(pii), copy(sigma), w,
           copy(tq))

    w = w[desc.ro, :]
    z = broadcast(.+, broadcast(.*, z, desc.sdv'), desc.xbar')
    z = z[desc.ro, :]

    zorig = data[:, (desc.p + 1):]
    z[zorig .!= mvcode] = zorig[zorig .!= mvcode]
    w, z
end

end # module
