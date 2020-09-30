"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored and optimized by Eric Ford
"""

default_measure_width_at_frac_depth = 0.5
default_frac_of_width_to_fit = 0.5
default_init_guess_ccf_σ = 2000.0   # m/s

# Abstract types
""" Abstract type for functors to estimate the raidal velocitiy from a CCF and its velocity grid.  """
abstract type AbstractMeasureRvFromCCF end

"""
Estimate RV based on centroid velocity of the CCF.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
"""
function (::AbstractMeasureRvFromCCF)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} } end

# Functions to be exported

"""
    `measure_rv_from_ccf(vels, ccf; alg )`
Return estimated RV based on the CCF using specified algorithm object.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Optional Arguements:
* `alg`: Functor specifying how to measure the RV from the CCF.  Options include:
MeasureRvFromCCFGaussian (default), MeasureRvFromCCFQuadratic, MeasureRvFromCCFCentroid, and MeasureRvFromMinCCF.
"""
function measure_rv_from_ccf(vels::A1, ccf::A2, ; alg::AlgT=MeasureRvFromCCFGaussian() )  where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, AlgT<:AbstractMeasureRvFromCCF }
    return alg(vels,ccf)
end

"""   `measure_rv_from_ccf(vels, ccf; alg )`
At each time, return the estimated radial velocity based on the CCF using the specified algorithm.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Optional Arguements:
* `alg`: Functor specifying how to measure the RV from the CCF.  Options include:
MeasureRvFromCCFGaussian (default), MeasureRvFromCCFQuadratic, MeasureRvFromCCFCentroid, and MeasureRvFromMinCCF.
"""
function measure_rvs_from_ccf(vels::A1, ccf::A2; alg::AlgT=MeasureRvFromCCFGaussian() )  where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, AlgT<:AbstractMeasureRvFromCCF }
    RVs = zeros(size(ccf,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv_from_ccf(vels, view(ccf,:,i), alg=alg)
    end
    return RVs
end

# Utility functions

"""  `est_full_width(vels, ccf; measure_width_at_frac_depth = 0.5 )`
Return rough estimate of a ccf (or line) full width at the specified fractional depth (i.e., fraction of total depth).
This is based on the velocities of the first/last pixel to have a value less than the target value, with no interpolation.
Assumes vels is sorted.  Depth is measured assuming to the maximum value of ccf represents the continuum.
Could be improved for noisy data and segments with a steep slope due to the blaze or neighboring lines.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Optional Arguements:
* `measure_width_at_frac_depth`: What fractional CCF depth should be used for defining the width (0.5)
"""
function est_full_width(vels::A1, ccf::A2; measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth ) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    minccf, maxccf = extrema(ccf)
    depth = maxccf - minccf
    target_val = minccf + (1-measure_width_at_frac_depth) * depth
    ind1 = findfirst(ccf .<= target_val)
    ind2 = findlast(ccf .<= target_val)
    if isnothing(ind1) || isnothing(ind2)
        println("ccf = ",ccf)
        println("minccf= ", minccf, " maxccf= ", maxccf, " depth= ", depth, " measure_width_at_frac_depth= ", measure_width_at_frac_depth, " targetval= ",target_val, " ind1= ", ind1, " ind2= ", ind2)
        @error "est_full_width failed."
    end
    return vels[ind2] - vels[ind1]
end

"""  `find_idx_at_and_around_minimum(vels, ccf; frac_of_width, measure_width_at_frac_depth )`
Return the a pair with the index of the lowest value of ccf and a range of pixels surrounding it.
The range is based on finding the pixels with velocities nearest to the
Assumes vels is sorted.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Optional Arguements:
* `frac_of_width_to_fit`: How large of a range of velocities should be included in the fit (0.5)
* `measure_width_at_frac_depth`: What fractional CCF depth should be used for defining the width (0.5)
"""
function find_idx_at_and_around_minimum(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit, measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # do a prelim fit to get the width
    full_width = est_full_width(vels, ccf, measure_width_at_frac_depth=measure_width_at_frac_depth)

    # find the min and fit only that
    amin = argmin(ccf)
    lend = vels[amin] - frac_of_width_to_fit * full_width
    rend = vels[amin] + frac_of_width_to_fit * full_width

    # get the indices
    lind = searchsortednearest(vels[1:amin], lend)
    rind = amin + searchsortednearest(vels[amin+1:end], rend)
    inds = lind:rind

    return (amin, inds)
end
