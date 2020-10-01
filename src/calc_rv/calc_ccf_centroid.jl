"""
    Code to estimated RV from centroid of core of CCF
Author: Eric Ford
Created: September 2020
"""

"""  Functor to estimate RV based on the centroid of the CCF.  """
struct MeasureRvFromCCFCentroid <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* `frac_of_width_to_fit`: (0.5)
* `measure_width_at_frac_depth`: (0.5)
"""
function MeasureRvFromCCFCentroid(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 2.0
    MeasureRvFromCCFCentroid(frac_of_width_to_fit,measure_width_at_frac_depth)
end

#=
"""
Estimate RV based on centroid velocity of the CCF.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Only makes sense if the CCF has been variance normalized.  So don't use this version without ccf variances.
"""
function (mrv::MeasureRvFromCCFCentroid)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and use only that part of the CCF for computing centroid
    amin, inds = find_idx_at_and_around_minimum(vels, ccf,
                        frac_of_width_to_fit=mrv.frac_of_width_to_fit,
                        measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

    weights = exp.(-0.5.*(ccf[inds].-ccf[amin]))
    v_centroid = sum(weights.*vels[inds]) / sum(weights)
    return (rv=v_centroid, σ_rv=NaN)
end
=#

"""
Estimate RV based on centroid velocity of the CCF.
Inputs:
* `vels`: Array of velocites where CCF was evaluated.
* `ccf`:  Array of values of CCF
Only makes sense if the CCF has been variance normalized.  So don't use this yet.
"""
function (mrv::MeasureRvFromCCFCentroid)(vels::A1, ccf::A2, ccf_var::A3 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1} }
    # find the min and use only that part of the CCF for computing centroid
    amin, inds = find_idx_at_and_around_minimum(vels, ccf,
                        frac_of_width_to_fit=mrv.frac_of_width_to_fit,
                        measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

    weights = exp.(-0.5.*(ccf[inds].-ccf[amin])) ./ ccf_var[inds]
    v_centroid = sum(weights.*vels[inds]) / sum(weights)
    return (rv=v_centroid, σ_rv=NaN)
end
