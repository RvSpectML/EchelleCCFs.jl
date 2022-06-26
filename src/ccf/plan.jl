"""
    Code to specifing plan for calculating a CCF
Author: Eric Ford
Created: August 2020
"""

"""A struct implementing a specific plans describing where the CCF is to be evaluated should be a subtype of AbstractCCFPlan. """
abstract type AbstractCCFPlan end

""" Basic plan for computing the CCF roughly between v_center-v_max and v_center+v_max with step size v_step. """
mutable struct BasicCCFPlan{MST<:AbstractCCFMaskShape, LLT<:AbstractLineList} <: AbstractCCFPlan
    v_center::Float64
    v_step::Float64
    v_max::Float64
    v_range_no_mask_change::Float64
    mask_shape::MST
    line_list::LLT
    allow_nans::Bool

    function BasicCCFPlan(midpoint::Real,step::Real, max::Real, range_no_mask_change::Real,
                          mask_shape::MST, line_list::LLT, allow_nans::Bool ) where { MST<:AbstractCCFMaskShape, LLT<:AbstractLineList }
        @assert 1.0e3 <= max <=100.0e3    # Reasonable range m/s, designed to prevent mistakes
        @assert 1 < step < 1000           # Reasonable range m/s
        #@assert abs(midpoint) < max      # Not true for SOAP simulations to ESPRESSO G2 mask
        new{MST,LLT}(midpoint, step, max, range_no_mask_change, mask_shape, line_list, allow_nans)
    end

end

"""   BasicCCFPlan
# Optional arguments:
- `midpoint`: (`default_v_center`)
- `step`: (`default_v_step`)
- `max`: (`default_v_max`)
"""
function BasicCCFPlan(;midpoint::Real=default_v_center, step::Real=default_v_step, max::Real=default_v_max,
                       range_no_mask_change::Real=max, mask_shape::MST,
                       line_list::LLT, allow_nans::Bool = true ) where { MST<:AbstractCCFMaskShape, LLT<:AbstractLineList }
    BasicCCFPlan(midpoint, step, max, range_no_mask_change, mask_shape, line_list, allow_nans)
end

""" `calc_ccf_v_grid( plan )`
Return range where CCF is to be evaluated,
Centered around plan.v_center going up to at least plan.v_max in steps of size plan.v_step.
Units based on those in plan.
"""
function calc_ccf_v_grid(p::PlanT where PlanT<:BasicCCFPlan )
    #n = ceil(Int, p.v_max/p.v_step)
    n = ceil(Int, (p.v_max-p.v_center)/p.v_step)
    range(p.v_center-n*p.v_step, p.v_center+n*p.v_step, length=2*n+1)
end

""" `calc_length_ccf_v_grid( plan )`
Return number of points in the velocity grid (without needing to create the range).
"""
function calc_length_ccf_v_grid(p::PlanT where PlanT<:BasicCCFPlan )
    #n = ceil(Int, p.v_max/p.v_step)
    n = ceil(Int, (p.v_max-p.v_center)/p.v_step)
    length=2*n+1
end

function get_mask_shape(p::PlanT ) where { PlanT<:BasicCCFPlan }
    return p.mask_shape
end

function set_mask_shape!(p::PlanT, m::MST ) where { PlanT<:BasicCCFPlan, MST<:AbstractCCFMaskShape }
    p.mask_shape = m
    return p
end

function increase_mask_fwhm!(p::PlanT, Δfwhm::Real ) where { PlanT<:BasicCCFPlan }
    if Δfwhm > 0
        set_mask_shape!(p, mask_with_increased_fwhm(get_mask_shape(p),Δfwhm) )
    end
    return p
end

function copy(plan::BasicCCFPlan{MST, LLT} ) where { MST<:AbstractCCFMaskShape, LLT<:AbstractLineList }
    BasicCCFPlan(plan.v_center, plan.v_step, plan.v_max, plan.v_range_no_mask_change, deepcopy(plan.mask_shape), plan.line_list, plan.allow_nans)
end
