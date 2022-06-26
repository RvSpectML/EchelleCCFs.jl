"""
    Code for convenience functions for calculating CCFs for an AbstractChunkList (see RvSpectMLBase)
Author: Eric Ford
Created: August 2020
"""

"""  `calc_ccf_chunklist ( chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                Δfwhm::Real = 0,
                                assume_sorted::Bool = false ) where {
                                            PlanT<:AbstractCCFPlan }
   num_chunks = length(chunk_list)
   @assert length(plan_for_chunk) == num_chunks
   num_vels = calc_length_ccf_v_grid(first(plan_for_chunk))
   ccf_out = Array{Float64,1}(undef, num_vels)
   if Δfwhm > 0
      this_plan_for_chunk = copy(plan_for_chunk)
      #increase_mask_fwhm!(this_plan_for_chunk,Δfwhm)
      map(chid->increase_mask_fwhm!(this_plan_for_chunk[chid],Δfwhm), 1:num_chunks )
   else
      this_plan_for_chunk = plan_for_chunk
   end
   calc_ccf_chunklist!(ccf_out, chunk_list, this_plan_for_chunk, assume_sorted=assume_sorted )
   return ccf_out
end

"""  `calc_ccf_chunklist! ( ccfs_out, chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist!(ccf_out::AbstractArray{T1,1}, chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                assume_sorted::Bool = false ) where {
                                            T1<:Real, PlanT<:AbstractCCFPlan }
   num_chunks = length(chunk_list)
   @assert length(plan_for_chunk) == num_chunks
   num_vels = calc_length_ccf_v_grid(first(plan_for_chunk))
   @assert size(ccf_out,1) == num_vels
   @assert all(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ) .== num_vels )
   @assert assume_sorted || all(map(chid->issorted( plan_for_chunk[chid].line_list.λ ), 1:num_chunks))

   num_lines = mapreduce(chid->length(plan_for_chunk[chid].line_list), +, 1:length(chunk_list.data) )
   ccf_out .= mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid],
                    assume_sorted=true ), +, 1:length(chunk_list.data) )
   ccf_out .*= num_lines
end

"""  `calc_ccf_and_var_chunklist ( chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_and_var_chunklist(chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                ccf_var_scale::Real = 1.0, Δfwhm::Real = 0,
                                assume_sorted::Bool = false ) where {
                                            PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  num_chunks = length(chunk_list)
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  ccf_out = Array{Float64,1}(undef, num_vels)
  ccf_var_out = Array{Float64,1}(undef, num_vels)

  if Δfwhm > 0
     this_plan_for_chunk = copy(plan_for_chunk)
     map(chid->increase_mask_fwhm!(this_plan_for_chunk[chid],Δfwhm), 1:num_chunks )
  else
     this_plan_for_chunk = plan_for_chunk
  end

  #(ccf_out, ccf_var_out ) = mapreduce(chid->calc_ccf_and_var_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                 assume_sorted=assume_sorted  ), add_tuple_sum, 1:length(chunk_list.data) )
  calc_ccf_and_var_chunklist!(ccf_out, ccf_var_out, chunk_list,this_plan_for_chunk, ccf_var_scale=ccf_var_scale, assume_sorted=assume_sorted)

  return (ccf=ccf_out, ccf_var=ccf_var_out)
end

"""  `calc_ccf_and_var_chunklist! ( ccf_out, ccf_var_out, chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_and_var_chunklist!(ccfs_out::AbstractArray{T1,1}, ccf_vars_out::AbstractArray{T2,1},
                                chunk_list::AbstractChunkList, plan_for_chunk::AbstractVector{PlanT};
                                ccf_var_scale::Real = 1.0, assume_sorted::Bool = false ) where {
                                            T1<:Real, T2<:Real, PlanT<:AbstractCCFPlan }
  num_chunks = length(chunk_list)
  @assert length(plan_for_chunk) == num_chunks
  @assert assume_sorted || all(map(chid->issorted( plan_for_chunk[chid].line_list.λ ), 1:num_chunks))
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  @assert size(ccfs_out,1) == num_vels
  @assert size(ccf_vars_out,1) == num_vels

  #(ccf_out, ccf_var_out ) = mapreduce(chid->calc_ccf_and_var_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                  assume_sorted=assume_sorted  ), add_tuple_sum, 1:length(chunk_list.data) )

  ccfs_out .= 0
  ccf_vars_out .= 0
  ccfs_tmp = zeros(num_vels)
  ccf_vars_tmp = zeros(num_vels)
  for chid in 1:num_chunks
     calc_ccf_and_var_chunk!(ccfs_tmp, ccf_vars_tmp,
        chunk_list.data[chid], plan_for_chunk[chid], ccf_var_scale=ccf_var_scale, assume_sorted=true )
     ccfs_out .+= ccfs_tmp
     ccf_vars_out .+= ccf_vars_tmp
  end

  return (ccfs=ccfs_out, ccf_vars=ccf_vars_out)
end

"""  `calc_ccf_and_covar_chunklist! ( ccf_out, ccf_cocovar_out, chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_and_covar_chunklist!(ccfs_out::AbstractArray{T1,1}, ccf_covars_out::AbstractArray{T2,2},
                                chunk_list::AbstractChunkList, plan_for_chunk::AbstractVector{PlanT};
                                ccf_var_scale::Real = 1.0, assume_sorted::Bool = false ) where {
                                            T1<:Real, T2<:Real, PlanT<:AbstractCCFPlan }
  num_chunks = length(chunk_list)
  @assert length(plan_for_chunk) == num_chunks
  @assert assume_sorted || all(map(chid->issorted( plan_for_chunk[chid].line_list.λ ), 1:num_chunks))
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  @assert size(ccfs_out,1) == num_vels
  @assert size(ccf_covars_out,1) == size(ccf_covars_out,2) == num_vels

  #(ccf_out, ccf_covar_out ) = mapreduce(chid->calc_ccf_and_covar_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                  assume_sorted=assume_sorted  ), add_tuple_sum, 1:length(chunk_list.data) )

  ccfs_out .= 0
  ccf_covars_out .= 0
  ccfs_tmp = zeros(num_vels)
  ccf_covars_tmp = zeros(num_vels,num_vels)
  for chid in 1:num_chunks
     calc_ccf_and_covar_chunk!(ccfs_tmp, ccf_covars_tmp,
        chunk_list.data[chid], plan_for_chunk[chid], ccf_var_scale=ccf_var_scale, assume_sorted=true )
     ccfs_out .+= ccfs_tmp
     ccf_covars_out .+= ccf_covars_tmp
  end

  return (ccfs=ccfs_out, ccf_covars=ccf_covars_out)
end

"""  `calc_ccf_and_covar_chunklist ( chunklist, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_and_covar_chunklist(chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                ccf_var_scale::Real = 1.0,
                                Δfwhm::Real = 0,
                                assume_sorted::Bool = false ) where {
                                            PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  num_chunks = length(chunk_list)
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  ccf_out = Array{Float64,1}(undef, num_vels)
  ccf_covar_out = Array{Float64,2}(undef, num_vels, num_vels)

  if Δfwhm > 0
     this_plan_for_chunk = copy(plan_for_chunk)
     #increase_mask_fwhm!(this_plan_for_chunk,Δfwhm)
     map(chid->increase_mask_fwhm!(this_plan_for_chunk[chid],Δfwhm), 1:num_chunks )
  else
     this_plan_for_chunk = plan_for_chunk
  end

  #(ccf_out, ccf_covar_out ) = mapreduce(chid->calc_ccf_and_covar_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                 assume_sorted=assume_sorted  ), add_tuple_sum, 1:length(chunk_list.data) )
  calc_ccf_and_covar_chunklist!(ccf_out, ccf_covar_out, chunk_list,this_plan_for_chunk, ccf_var_scale=ccf_var_scale, assume_sorted=assume_sorted)

  return (ccf=ccf_out, ccf_covar=ccf_covar_out)
end
