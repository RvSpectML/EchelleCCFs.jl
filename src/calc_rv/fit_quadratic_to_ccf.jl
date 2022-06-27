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

using Statistics
function (mrv::MeasureRvFromCCFQuadratic)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
    if isnan(amin)     return ( rv=NaN, σ_rv=NaN )     end

    mean_v = mean(view(vels,inds))
    X = ones(length(inds),3)
    X[:,2] .= view(vels,inds) .-mean_v
    X[:,3] .= (view(vels,inds) .-mean_v).^2
    (c, b, a)  = (X'*X) \ (X'*view(ccf,inds))

    v_at_min_of_quadratic = -b/(2*a) + mean_v
    return (rv=v_at_min_of_quadratic, σ_rv=NaN)
end

function (mrv::MeasureRvFromCCFQuadratic)(vels::A1, ccf::A2, ccf_var::A3 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1} }
    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
    if isnan(amin)     return ( rv=NaN, σ_rv=NaN )     end
    mean_v = mean(view(vels,inds))
    X = ones(length(inds),3)
    X[:,2] .= view(vels,inds) .-mean_v
    X[:,3] .= (view(vels,inds) .-mean_v).^2
    covar = PDiagMat(abs.(view(ccf_var,inds)))
    denom = (X' * (covar \ X) )
    (c, b, a)  = denom \ (X' * (covar \ view(ccf,inds)) )
    covar_beta = inv(denom)
    v_at_min_of_quadratic = -b/(2*a)
    sigma_rv = abs(v_at_min_of_quadratic) * sqrt(covar_beta[2,2]/b^2+covar_beta[3,3]/a^2)
    v_at_min_of_quadratic += mean_v
    return (rv=v_at_min_of_quadratic, σ_rv=sigma_rv)
end
