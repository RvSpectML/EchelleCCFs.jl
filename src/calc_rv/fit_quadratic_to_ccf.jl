"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactors & optimized by Eric Ford
"""

"""  Functor to estimate RV based on fitting quadratic near minimum of CCF.
TODO: Revist the logic here and see if need to perform transformation first.
"""
struct MeasureRvFromCCFQuadratic <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* `frac_of_width_to_fit`: (0.5)
* `measure_width_at_frac_depth`: (0.5)
"""
function MeasureRvFromCCFQuadratic(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 2.0
    MeasureRvFromCCFQuadratic(frac_of_width_to_fit,measure_width_at_frac_depth)
end

function (mrv::MeasureRvFromCCFQuadratic)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

    # do the polyfit
    pfit = Polynomials.fit(vels[inds], ccf[inds], 2)
    @assert length(Polynomials.coeffs(pfit)) >= 3   # just in case fails to fit a quadratic

    # get center from coeffs
    c, b, a = Polynomials.coeffs(pfit)
    v_at_min_of_quadratic = -b/(2*a)
    return v_at_min_of_quadratic
end
