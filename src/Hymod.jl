__precompile__()

module Hymod

using DataFrames
using Dates
using Distributions
using Statistics
using Logging
using JSON

export randomparams, hargreaves, simulate, calibrate

"""
    randomparams()

Function to create a random set of parameters as kwargs for Hymod.simulate()
"""
function randomparams()

    param_bounds = Dict{Symbol,Dict}(
        :cmax => Dict{Symbol,Float64}(:lower => 1.0, :upper => 100),
        :bexp => Dict{Symbol,Float64}(:lower => 0.0, :upper => 2.0),
        :alpha => Dict{Symbol,Float64}(:lower => 0.2, :upper => 0.99),
        :ks => Dict{Symbol,Float64}(:lower => 0.01, :upper => 0.5),
        :kq => Dict{Symbol,Float64}(:lower => 0.5, :upper => 1.2),
    )

    out = Dict{Symbol,Float64}()

    for k in keys(param_bounds)
        minV = param_bounds[k][:lower]
        maxV = param_bounds[k][:upper]
        p = rand(Uniform(minV, maxV), 1)
        out[k] = p[1]
    end

    return out
end

function abspower(x::Float64, y::Float64)
    # Needed to capture invalid overflow with netgative values
    x = abs(x)
    return x^y
end

function precipexcess(
    x_loss::Float64,
    cmax::Float64,
    bexp::Float64,
    prc::Float64,
    pet::Float64,
)
    xn_prev = x_loss
    ct_prev =
        cmax *
        (1 - abspower((1 - ((bexp + 1) * (xn_prev) / cmax)), (1 / (bexp + 1))))
    # Calculate Effective rainfall 1
    ER1 = max((prc - cmax + ct_prev), 0.0)
    prc = prc - ER1
    dummy = min(((ct_prev + prc) / cmax), 1)
    xn = (cmax / (bexp + 1)) * (1 - abspower((1 - dummy), (bexp + 1)))

    # Calculate Effective rainfall 2
    ER2 = max(prc - (xn - xn_prev), 0)

    # Alternative approach
    # actual ET is linearly related to the soil moisture state
    evap = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * pet
    xn = max(xn - evap, 0)  # update state

    return ER1, ER2, xn
end


"""
    randomparams(x_slow::Float64, inflow::Float64, Rs::Float64)

Linear reservoir routing function
"""
function linearreservoir(x_slow::Float64, inflow::Float64, Rs::Float64)
    # Linear reservoir routing
    x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
    outflow = (Rs / (1 - Rs)) * x_slow
    return x_slow, outflow
end


"""
    hargreaves(forcings::DataFrame; tminCol::Symbol = :tmin, tmaxCol::Symbol = :tmax, dtCol::Symbol = :datetime)

Function to calculate PET using Hargreves equation.
"""
function hargreaves(
    forcings::DataFrame;
    tminCol::Symbol = :tmin,
    tmaxCol::Symbol = :tmax,
    dtCol::Symbol = :datetime,
)

    dts = forcings[!, dtCol]
    tmin = forcings[!, tminCol]
    tmax = forcings[!, tmaxCol]

    hargreaves(tmin, tmax, dts)
end

"""
    hargreaves(tmin::AbstractArray, tmax::AbstractArray, dts::AbstractArray)

Function to calculate PET using Hargreves equation.
"""
function hargreaves(
    tmin::AbstractArray,
    tmax::AbstractArray,
    dts::AbstractArray,
)
    hargreaves.(tmin, tmax, dts)
end


"""
    hargreaves(tmin::Real, tmax::Real, dts::Date)

Function to calculate PET using Hargreves equation.
"""
function hargreaves(tmin::Real, tmax::Real, dts::Date)
    Gsc = 367
    lhov = 2.257

    doy = dayofyear(dts)

    tavg = mean((tmin, tmax))

    b = 2 * pi * (Int16(doy) / 365)
    Rav =
        1.00011 +
        0.034221 * cos(b) +
        0.00128 * sin(b) +
        0.000719 * cos(2 * b) +
        0.000077 * sin(2 * b)
    Ho = ((Gsc * Rav) * 86400) / 1e6

    eto = (0.0023 * Ho * (tmax - tmin)^0.5 * (tavg + 17.8))

    return eto
end


"""
    simulate(forcings::DataFrame; precipCol::Symbol = :precip, petCol::Symbol = :pet, kwargs...)

Function to simulate discharge using Hymod rainfall-runoff model.
"""
function simulate(
    forcings::DataFrame;
    precipCol::Symbol = :precip,
    petCol::Symbol = :pet,
    kwargs...,
)
    p = forcings[:, precipCol]
    e = forcings[:, petCol]

    simulate(p, e; kwargs...)

end


"""
    simulate(precip::AbstractArray, pet::AbstractArray; initFlow::Bool = true, cmax::Float64 = 1.0, bexp::Float64 = 0.0, alpha::Float64 = 0.2, ks::Float64 = 0.01, kq::Float64 = 0.5)

Function to simulate discharge using Hymod rainfall-runoff model
"""
function simulate(
    precip::AbstractArray,
    pet::AbstractArray;
    initFlow::Bool = true,
    cmax::Float64 = 1.0,
    bexp::Float64 = 0.0,
    alpha::Float64 = 0.2,
    ks::Float64 = 0.01,
    kq::Float64 = 0.5,
)

    lt_to_m::Float64 = 0.001

    n_iters = length(precip)

    # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
    x_loss::Float64 = 0.0
    # Initialize slow tank state
    # value of 0 init flow works ok if calibration data starts with low discharge
    x_slow::Float64 = initFlow ? (2.3503 / (ks * 22.5)) : 0.0
    # Initialize state(s) of quick tank(s)
    x_quick = zeros(Float64, 3)
    outflow = zeros(Float64, n_iters)
    output = zeros(Float64, n_iters)
    # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

    # prealloc values within loops
    Pval::Float64 = 0.0
    PETval::Float64 = 0.0
    ER1::Float64 = 0.0
    ER2::Float64 = 0.0
    UQ::Float64 = 0.0
    US::Float64 = 0.0
    QS::Float64 = 0.0
    inflow::Float64 = 0.0

    @inbounds for i = 1:n_iters
        Pval = precip[i]
        PETval = pet[i]
        # Compute excess precipitation and evaporation
        ER1, ER2, x_loss = precipexcess(x_loss, cmax, bexp, Pval, PETval)
        # Calculate total effective rainfall
        ET = ER1 + ER2
        #  Now partition ER between quick and slow flow reservoirs
        UQ = alpha * ET
        US = (1 - alpha) * ET
        # Route slow flow component with single linear reservoir
        x_slow, QS = linearreservoir(x_slow, US, ks)
        # Route quick flow component with linear reservoirs
        inflow = UQ

        for j in eachindex(x_quick)
            # Linear reservoir
            x_quick[j], outflow = linearreservoir(x_quick[j], inflow, kq)
            inflow = outflow
        end

        # Compute total flow for timestep
        output[i] = ((QS + outflow) / lt_to_m)
    end

    output[output.<0] .= 0

    return output
end


"""
    calibrate(forcing::DataFrame, paramSpace::Dict, nSamples::Int64; precipCol::Symbol = :precip, petCol::Symbol = :pet, obsCol::Symbol = :obs, savefinal::Bool = false)

Function to calibrate Hymod rainfall-runoff model
"""
function calibrate(
    forcing::DataFrame,
    paramSpace::Dict,
    nSamples::Int64;
    precipCol::Symbol = :precip,
    petCol::Symbol = :pet,
    obsCol::Symbol = :obs,
    savefinal::Bool = false,
)
    # implementation of a monte carlo sampling

    obs = forcing[!, obsCol]

    params = Dict{Symbol,Array}()

    keyList = collect(keys(paramSpace))

    for k in keyList
        p = paramSpace[k]
        minV = p[:lower]
        maxV = p[:upper]
        params[k] = collect(rand(Uniform(minV, maxV), nSamples))
    end

    nIter = nSamples
    losses = zeros(nIter)

    @info "Calibrating model with $nIter iterations..."
    for i = 1:nIter
        vals = [params[k][i] for k in keyList]
        pars = Dict(keyList .=> vals)
        q = simulate(forcing, precipCol = :precip, petCol = :pet; pars...)
        # only test after one year of spinup
        losses[i] = nse(q[365:end], obs[365:end])
    end

    # check to make sure that the loss metric is a valid number
    losses = convert(Vector{Union{Missing,Float64}}, losses)
    losses[map(isnan, losses)] .= missing

    finalLoss, idx = findmax(skipmissing(losses))
    vals = [params[k][idx] for k in keyList]
    finalPars = Dict(keyList .=> vals)
    finalQ = simulate(forcing, precipCol = :precip, petCol = :pet; finalPars...)

    if savefinal
        t = now()
        finalPars[:loss] = finalLoss
        jsondata = Dict(String(k) => v for (k, v) in finalPars)
        @debug (jsondata)
        open("hymod_calibration_results_$t.json", "w") do io
            JSON.print(io, jsondata, 4)
        end
    end

    return finalQ, finalPars, finalLoss

end

"""
    nse(y_true::AbstractArray, y_pred::AbstractArray)

Function to calculate Nash-Sutcliffe model efficiency coefficient
"""
function nse(y_true::AbstractArray, y_pred::AbstractArray)
    # sum square difference between observed and simulated
    numerator = sum((y_true - y_pred) .^ 2)
    # sum square difference between observed and observed mean
    denominator = sum((y_true .- mean(y_true)) .^ 2)

    return 1 - (numerator / denominator)
end


end # module
