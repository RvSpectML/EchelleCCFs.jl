"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactors & optimized by Eric Ford
"""

#module FitGaussianToCCF

#import ..RVFromCCF: AbstractMeasureRvFromCCF, MeasureRvFromCCFQuadratic
# packages
using LsqFit
using LinearAlgebra

#export MeasureRvFromCCFGaussian

"""  Functor to estimate RV based on fitting a Gaussian quadratic near minimum of the CCF. """
struct MeasureRvFromCCFGaussian <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
    init_guess_ccf_σ::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* `frac_of_width_to_fit`: (0.5)
* `measure_width_at_frac_depth`: (0.5)
* `init_guess_ccf_σ``: (2000m/s)
"""
function MeasureRvFromCCFGaussian(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth,
                                    init_guess_ccf_σ::Real = default_init_guess_ccf_σ )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 5.0  # TODO: FIgure out appropriate range
    @assert 1 <= init_guess_ccf_σ <= 30e3
    MeasureRvFromCCFGaussian(frac_of_width_to_fit,measure_width_at_frac_depth,init_guess_ccf_σ)
end

@. gaussian_line_helper(x, p) = p[4] + p[3] * exp(-0.5*((x-p[1])/p[2])^2)

function (mrv::MeasureRvFromCCFGaussian)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
        # find the min and fit only the part near the minimum of the CCF
        amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

        # make initial guess parameters
        μ = vels[amin]
        σ = mrv.init_guess_ccf_σ                          # TODO: Make sure this is robust.
        minccf, maxccf = extrema(ccf)
        amp = minccf - maxccf
        y0 = maxccf
        p0 = [μ, σ, amp, y0]

        # fit and return the mean of the distribution
        result = curve_fit(gaussian_line_helper, view(vels,inds), view(ccf,inds), p0)

        if result.converged
           rv = coef(result)[1]
           sigma_rv = stderror(result)[1]
           rvfit = (rv=rv, σ_rv=sigma_rv)
        else
           @warn "Fit of Gaussian to CCF Failed.  Reverting to fit quadratic to CCF."
           quad_fit_to_ccf = MeasureRvFromCCFQuadratic(frac_of_width_to_fit=frac_of_width_to_fit,measure_width_at_frac_depth=measure_width_at_frac_depth)
           rvfit = quad_fit_to_ccf(vels,ccf)
        end
        return rvfit
end


function (mrv::MeasureRvFromCCFGaussian)(vels::A1, ccf::A2, ccf_var::A3 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1} }
        # find the min and fit only the part near the minimum of the CCF
        amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

        # make initial guess parameters
        μ = vels[amin]
        σ = mrv.init_guess_ccf_σ                          # TODO: Make sure this is robust.
        minccf, maxccf = extrema(ccf)
        amp = minccf - maxccf
        y0 = maxccf
        p0 = [μ, σ, amp, y0]

        # fit and return the mean of the distribution
        result = curve_fit(gaussian_line_helper, view(vels,inds), view(ccf,inds), 1.0 ./ view(ccf_var,inds),  p0)

        if result.converged
           rv = coef(result)[1]
           sigma_rv = stderror(result)[1]
           rvfit = (rv=rv, σ_rv=sigma_rv)
        else
           @warn "Fit of Gaussian to CCF Failed.  Reverting to fit quadratic to CCF."
           quad_fit_to_ccf = MeasureRvFromCCFQuadratic(frac_of_width_to_fit=mrv.frac_of_width_to_fit,measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
           rvfit = quad_fit_to_ccf(vels,ccf,ccf_var)
        end
        return rvfit
end

#end # module FitGaussianToCCF
