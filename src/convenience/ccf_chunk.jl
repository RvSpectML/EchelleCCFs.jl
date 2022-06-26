"""
    Code for convenience functions for calculating CCFs for an AbstractChunkOfSpectrum (see RvSpectMLBase)
Author: Eric Ford
Created: August 2020
"""

"""  `calc_ccf_chunk!( ccf_out, chunk, ccf_plan )`
Convenience function to compute CCF for one chunk of spectrum, evaluated using mask_shape and line list from ccf plan

# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Returns:
- `ccf_out`
"""
function calc_ccf_chunk!(ccf_out::AbstractArray{T1,1}, chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T2} = chunk.var,
                 assume_sorted::Bool = false ) where { T1<:Real, T2<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  @assert length(ccf_out) == calc_length_ccf_v_grid(plan)
  ccf_1D!(ccf_out, chunk.λ, chunk.flux, plan; assume_sorted=true)
  return ccf_out
end

"""  `calc_ccf_and_var_chunk!( chunk, ccf_plan )`
Convenience function to compute CCF and variance of each "CCF pixel" for one chunk of spectrum, evaluated using mask_shape and line list from `ccf_plan`.
# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `ccf_var_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `var`:  `AbstractArray` with variance to use for each pixel (overides value in chunk)
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Returns Named Tuple with:
- `ccf_out`:
- `ccf_var_out`:
"""
function calc_ccf_and_var_chunk!(ccf_out::AbstractArray{T1,1}, ccf_var_out::AbstractArray{T2,1},
                chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T3} = chunk.var, ccf_var_scale::Real = 1.0,
                 assume_sorted::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  @assert length(ccf_out) == calc_length_ccf_v_grid(plan)
  ccf_1D!(ccf_out, ccf_var_out, chunk.λ, chunk.flux, var, plan; ccf_var_scale=ccf_var_scale, assume_sorted=true)
  return (ccf=ccf_out, ccf_var=ccf_var_out)
end

"""  `calc_ccf_and_covar_chunk!( chunk, ccf_plan )`
Convenience function to compute CCF and covariance of each pair of "CCF pixels" for one chunk of spectrum, evaluated using mask_shape and line list from `ccf_plan`.
# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `ccf_covar_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `var`:  `AbstractArray` with variance to use for each pixel (overides value in chunk)
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Returns Named Tuple with:
- `ccf_out`:
- `ccf_covar_out`:
"""
function calc_ccf_and_covar_chunk!(ccf_out::AbstractArray{T1,1}, ccf_covar_out::AbstractArray{T2,2},
                chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T3} = chunk.var, ccf_var_scale::Real = 1.0,
                 assume_sorted::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  @assert length(ccf_out) == calc_length_ccf_v_grid(plan)
  ccf_1D!(ccf_out, ccf_covar_out, chunk.λ, chunk.flux, var, plan; ccf_var_scale=ccf_var_scale, assume_sorted=true)
  return (ccf=ccf_out, ccf_covar=ccf_covar_out)
end


"""  `calc_ccf_chunk( chunk, ccf_plan )`
Convenience function to compute CCF for one chunk of spectrum.
# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
- `calc_ccf_var`:  if true also computes estimate of variance for each value of ccf
# Returns:
- CCF for one chunk of spectrum, evaluated using mask_shape and line list from ccf plan
"""
function calc_ccf_chunk(chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T} = chunk.var, ccf_var_scale::Real = 1.0, Δfwhm::Real = 0,
                 assume_sorted::Bool = false, calc_ccf_var::Bool = false ) where { T<:Real, PlanT<:AbstractCCFPlan }
  if Δfwhm > 0
     this_plan_for_chunk = copy(plan)
     increase_mask_fwhm!(this_plan_for_chunk,Δfwhm)
  else
     this_plan_for_chunk = plan
  end
  len_v_grid = calc_length_ccf_v_grid(plan)
  ccf_out = zeros(len_v_grid)
  if calc_ccf_var
    ccf_var_out = zeros(len_v_grid)
    return calc_ccf_and_var_chunk!(ccf_out, ccf_var_out, chunk, this_plan_for_chunk, var=var, ccf_var_scale=ccf_var_scale, assume_sorted=assume_sorted )
  else
    return calc_ccf_chunk!(ccf_out, chunk, this_plan_for_chunk, var=var, assume_sorted=assume_sorted )
  end
end

"""  `calc_ccf_and_var_chunk( chunk, ccf_plan )`
Convenience function to compute CCF and variance of each "CCF pixel" for one chunk of spectrum, evaluated using mask_shape and line list from `ccf_plan`.
# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `ccf_var_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `var`:  `AbstractArray` with variance to use for each pixel (overides value in chunk)
`- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Returns Named Tuple with:
- `ccf_out`:
- `ccf_var_out`:
"""
function calc_ccf_and_var_chunk(chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T} = chunk.var, ccf_var_scale::Real = 1.0, Δfwhm::Real = 0,
                 assume_sorted::Bool = false ) where { T<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  if Δfwhm > 0
     this_plan_for_chunk = copy(plan)
     increase_mask_fwhm!(this_plan_for_chunk,Δfwhm)
  else
     this_plan_for_chunk = plan
  end
  len_v_grid = calc_length_ccf_v_grid(plan)
  ccf_out = zeros(len_v_grid)
  ccf_var_out = zeros(len_v_grid)
  calc_ccf_and_var_chunk!(ccf_out, ccf_var_out, chunk, this_plan_for_chunk, var=var, ccf_var_scale=ccf_var_scale, assume_sorted=true )
  return (ccf=ccf_out, ccf_var=ccf_var_out)
end

"""  `calc_ccf_and_covar_chunk( chunk, ccf_plan )`
Convenience function to compute CCF and covariance of each pair of "CCF pixels" for one chunk of spectrum, evaluated using mask_shape and line list from `ccf_plan`.
# Inputs:
- `ccf_out`:  `AbstractArray` to store output
- `ccf_covar_out`:  `AbstractArray` to store output
- `chunk`: ChunkOfSpectrum to compute CCF for
- `ccf_plan`: for now, just a BasicCCFPlan that provides line_list, mask_shape and other parameters for calculating CCF
# Optional Arguments:
- `var`:  `AbstractArray` with variance to use for each pixel (overides value in chunk)
`- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Returns Named Tuple with:
- `ccf_out`:
- `ccf_covar_out`:
"""
function calc_ccf_and_covar_chunk(chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T} = chunk.var, ccf_var_scale::Real = 1.0, Δfwhm::Real = 0,
                 assume_sorted::Bool = false ) where { T<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  if Δfwhm > 0
     this_plan_for_chunk = copy(plan)
     increase_mask_fwhm!(this_plan_for_chunk,Δfwhm)
  else
     this_plan_for_chunk = plan
  end
  len_v_grid = calc_length_ccf_v_grid(plan)
  ccf_out = zeros(len_v_grid)
  ccf_covar_out = zeros(len_v_grid,len_v_grid)
  calc_ccf_and_covar_chunk!(ccf_out, ccf_covar_out, chunk, this_plan_for_chunk, var=var, ccf_var_scale=ccf_var_scale, assume_sorted=true )
  return (ccf=ccf_out, ccf_covar=ccf_covar_out)
end
