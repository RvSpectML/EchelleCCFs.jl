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
                                assume_sorted::Bool = false #=, use_pixel_vars::Bool = false =#  ) where {
                                            PlanT<:AbstractCCFPlan }
   num_vels = calc_length_ccf_v_grid(first(plan_for_chunk))
   ccf_out = Array{Float64,1}(undef, num_vels)
   calc_ccf_chunklist!(ccf_out, chunk_list, plan_for_chunk, assume_sorted=assume_sorted )
   return ccf_out
   #mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid],
   #                 assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =#), +, 1:length(chunk_list.data) )
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
                                assume_sorted::Bool = false #=, use_pixel_vars::Bool = false =#  ) where {
                                            T1<:Real, PlanT<:AbstractCCFPlan }
   num_chunks = length(chunk_list)
   @assert length(plan_for_chunk) == num_chunks
   num_vels = calc_length_ccf_v_grid(first(plan_for_chunk))
   @assert size(ccf_out,1) == num_vels
   @assert all(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ) .== num_vels )
   @assert assume_sorted || all(map(chid->issorted( plan_for_chunk[chid].line_list.λ ), 1:num_chunks))

   ccf_out .= mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid],
                    assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =#), +, 1:length(chunk_list.data) )
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
                                assume_sorted::Bool = false #=, use_pixel_vars::Bool = false =#  ) where {
                                            PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  num_chunks = length(chunk_list)
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  ccf_out = Array{Float64,1}(undef, num_vels)
  ccf_var_out = Array{Float64,1}(undef, num_vels)

  #(ccf_out, ccf_var_out ) = mapreduce(chid->calc_ccf_and_var_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                 assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =# ), add_tuple_sum, 1:length(chunk_list.data) )
  calc_ccf_and_var_chunklist!(ccf_out, ccf_var_out, chunk_list,plan_for_chunk, assume_sorted=assume_sorted)

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
                                assume_sorted::Bool = false #=, use_pixel_vars::Bool = false =#  ) where {
                                            T1<:Real, T2<:Real, PlanT<:AbstractCCFPlan }
  num_chunks = length(chunk_list)
  @assert length(plan_for_chunk) == num_chunks
  @assert assume_sorted || all(map(chid->issorted( plan_for_chunk[chid].line_list.λ ), 1:num_chunks))
  num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
  @assert size(ccfs_out,1) == num_vels

  #(ccf_out, ccf_var_out ) = mapreduce(chid->calc_ccf_and_var_chunk(chunk_list.data[chid], plan_for_chunk[chid],
  #                  assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =# ), add_tuple_sum, 1:length(chunk_list.data) )

  ccfs_out .= 0
  ccf_vars_out .= 0
  ccfs_tmp = zeros(num_vels)
  ccf_vars_tmp = zeros(num_vels)
  for chid in 1:num_chunks
     calc_ccf_and_var_chunk!(ccfs_tmp, ccf_vars_tmp,
        chunk_list.data[chid], plan_for_chunk[chid], assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =# )
     ccfs_out .+= ccfs_tmp
     ccf_vars_out .+= ccf_vars_tmp
  end

  return (ccfs=ccfs_out, ccf_vars=ccf_vars_out)
end


#=
"""  `calc_ccf_chunklist ( chunklist, var_list, ccf_plans )`
Convenience function to compute CCF based on a spectrum's chunklist.
Prototype/Experimental version trying to use pixel variances not yet fully implemented/tested.
# Inputs:
- chunklist
- var_list: vector of variance vectors for each chunk
- ccf_palns: vector of ccf plans (one for each chunk)
# Optional Arguments:
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList, var_list::AbstractVector{A1},
                                plan_for_chunk::AbstractVector{PlanT};
                                assume_sorted::Bool = false, use_pixel_vars::Bool = false   ) where {
                                            T1<:Real, A1<:AbstractVector{T1}, PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  @assert length(chunk_list) == length(var_list)
  @assert all(map(chid->length(chunk_list.data[chid].var) == length(var_list[chid]),1:length(chunk_list.data)))
  #var_weight = mapreduce(chid->median(chunk_list.data[chid].var), +, 1:length(chunk_list.data) ) / mapreduce(chid->median(var_list[chid]), +, 1:length(chunk_list.data) )
  var_weight = 1.0 #mapreduce(chid->mean(chunk_list.data[chid].var), +, 1:length(chunk_list.data) ) / mapreduce(chid->mean(var_list[chid]), +, 1:length(chunk_list.data) )
  #println("# var_weight*N_obs = ", var_weight .* length(chunk_list))
  mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid], var=var_list[chid] .* var_weight,
                    assume_sorted=assume_sorted, #= use_pixel_vars=use_pixel_vars), =# +, 1:length(chunk_list.data) )
end
=#
