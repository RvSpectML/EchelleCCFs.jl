"""
Author: Eric Ford
Created: September 2020
"""

#module FitTemplateToCCF

#import ..RVFromCCF: AbstractMeasureRvFromCCF, MeasureRvFromCCFQuadratic
# packages
using LinearAlgebra

# TODO: export once passes test that it's working
#export MeasureRvFromCCFTemplate

"""  Functor to estimate RV based on fitting a template and its derivative near minimum of the CCF. """
struct MeasureRvFromCCFTemplate <: AbstractMeasureRvFromCCF
    v_grid::AbstractVector{Float64}
    template::AbstractVector{Float64}
    deriv::AbstractVector{Float64}
    v_idx_to_fit::UnitRange{Int64}
end

# Todo: Replace with version in RvSpectMLBase or a brute force GP?
function numerical_deriv( x::AA1, y::AA2 )   where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}  }
	@assert length(x) == length(y)
	dfluxdlnλ = zeros(size(y))
	dfluxdlnλ[1] = (y[2]-y[1])/(x[2]-x[1])
	dfluxdlnλ[2:end-1] .= (y[3:end].-y[1:end-2])./(x[3:end].-x[1:end-2])
	dfluxdlnλ[end] = (y[end]-y[end-1])/(x[end]-x[end-1])
	return dfluxdlnλ
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
- `v_grid`:  list of velocities where template is evaluated
- `tempalte`: template ccf evaluated at v_grid
- `frac_of_width_to_fit`: (0.5)
- `measure_width_at_frac_depth`: (0.5)
"""
function MeasureRvFromCCFTemplate(; v_grid::AbstractVector{T1},
                                    template::AbstractVector{T2},
                                    frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth
                                     ) where { T1<:Real, T2<:Real }
    @assert length(v_grid) == length(template)
    @assert length(v_grid >= 3)
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 5.0  # TODO: FIgure out appropriate range
	# Estimate derivative
    deriv = numerical_deriv(v_grid, template)
	# find the min and fit only the part near the minimum of the CCF
    v_min, v_idx_to_fit = find_idx_at_and_around_minimum(v_grid, template, frac_of_width_to_fit=frac_of_width_to_fit, measure_width_at_frac_depth=measure_width_at_frac_depth)
    MeasureRvFromCCFTemplate(v_grid, template, deriv, v_idx_to_fit)
end

function (mrv::MeasureRvFromCCFTemplate)(vels::A1, ccf::A2; ccf_var::A3 = zeros(length(ccf)) ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1}  }
        @assert all(vels .== mrv.v_grid)  # Could relax this, but keep it simple for now
        # find the min and fit only the part near the minimum of the CCF
        norm = sum(abs2.(mrv.deriv[v_idx_to_fit]))
        rv = sum((ccf[v_idx_to_fit].-mrv.template[v_idx_to_fit]).*mrv.deriv[v_idx_to_fit],dims=1).*(speed_of_light_mps/norm)
        # TODO: WARN: Uncertinaties ignores correlations between pixels, particularly problematic when oversample pixels
        σ_rv = sqrt.(sum(ccf_var[v_idx_to_fit].*abs2.(deriv[v_idx_to_fit]))) .*(speed_of_light_mps/norm)
		return (rv=rv, σ_rv=σ_rv)
end
