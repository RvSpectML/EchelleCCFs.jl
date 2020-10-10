"""
    Code for convenience functions for calculating seprate CCFs for each chunk in an AbstractChunkList (see RvSpectMLBase)
Author: Eric Ford
Created: August 2020
"""

"""  `calc_order_ccfs_chunklist ( chunklist_timeseries, list_of_ccf_plans )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) in a spectrum.
CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
# Inputs:
- `chunklist_timeseries`:
- `list_of_ccf_plans`: ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
- `use_pixel_vars`:  Don't set this true (yet)
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccfs_chunklist!(chunk_ccfs_out::AbstractArray{T1,2}, chunk_list::AbstractChunkList,
    plan_for_chunk::AbstractVector{PlanT} = BasicCCFPlan(); assume_sorted::Bool = false, use_pixel_vars::Bool = false  ) where {
                        T1<:Real, PlanT<:AbstractCCFPlan }
    num_chunks = length(chunk_list)
    @assert length(plan_for_chunk) == num_chunks
    @assert size(chunk_ccfs_out,2) == num_chunks
    @assert assume_sorted || issorted( first(plan_for_chunk).line_list.λ )
    num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
    @assert size(chunk_ccfs_out,1) == num_vels
    #flush(stdout); println("num_vels = ",num_vels," num_chunks = ", num_chunks, " length(chunk_list.data)=", length(chunk_list.data));     flush(stdout);
    for chid in 1:num_chunks
        # chunk_ccfs_out[:,i] .=
        #flush(stdout);         println("chid = ", chid);         flush(stdout);
        calc_ccf_chunk!(view(chunk_ccfs_out, 1:num_vels, chid ), chunk_list.data[chid], plan_for_chunk[chid], assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars)
    end
    return chunk_ccfs_out
end

"""  `calc_order_ccfs_chunklist ( chunklist_timeseries, list_of_ccf_plans )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) in a spectrum.
CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
# Inputs:
- `chunklist_timeseries`:
- `list_of_ccf_plans`: ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
- `use_pixel_vars`:  Don't set this true (yet)
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccfs_chunklist(chunk_list::AbstractChunkList,
    plan_for_chunk::AbstractVector{PlanT} = BasicCCFPlan(); assume_sorted::Bool = false, use_pixel_vars::Bool = false  ) where {
                        PlanT<:AbstractCCFPlan }
    num_chunks = length(chunk_list)
    @assert length(plan_for_chunk) == num_chunks
    @assert assume_sorted || issorted( first(plan_for_chunk).line_list.λ )
    num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
    chunk_ccfs_out = zeros( num_vels, num_chunks )
    calc_order_ccfs_chunklist!(chunk_ccfs_out, chunk_list, plan_for_chunk, assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars)
    return chunk_ccfs_out
end


"""  `calc_order_ccfs_chunklist ( chunklist_timeseries, list_of_ccf_plans )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) in a spectrum.
CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
# Inputs:
- `chunklist_timeseries`:
- `list_of_ccf_plans`: ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
- `use_pixel_vars`:  Don't set this true (yet)
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccf_and_vars_chunklist!(chunk_ccfs_out::AbstractArray{T1,2}, chunk_ccf_vars_out::AbstractArray{T1,2}, chunk_list::AbstractChunkList,
    plan_for_chunk::AbstractVector{PlanT} = BasicCCFPlan(); assume_sorted::Bool = false, use_pixel_vars::Bool = false  ) where {
                        T1<:Real, PlanT<:AbstractCCFPlan }
    num_chunks = length(chunk_list)
    @assert length(plan_for_chunk) == num_chunks
    @assert size(chunk_ccfs_out,2) == num_chunks
    @assert assume_sorted || issorted( first(plan_for_chunk).line_list.λ )
    num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
    @assert size(chunk_ccfs_out,1) == num_vels

    for chid in 1:num_chunks
        calc_ccf_and_var_chunk!(view(chunk_ccfs_out, 1:num_vels, chid ), view(chunk_ccf_vars_out, 1:num_vels, chid ),
            chunk_list.data[chid], plan_for_chunk[chid], assume_sorted=assume_sorted #=, use_pixel_vars=use_pixel_vars =# )
    end
    return (ccfs=chunk_ccfs_out, ccf_vars=chunk_ccf_vars_out)
end

"""  `calc_order_ccfs_chunklist ( chunklist_timeseries, list_of_ccf_plans )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) in a spectrum.
CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
# Inputs:
- `chunklist_timeseries`:
- `list_of_ccf_plans`: ccf plans (one for each chunk)
# Optional Arguments:
- `assume_sorted`:  if true, skips checking the line_list is sorted by wavelength
- `use_pixel_vars`:  Don't set this true (yet)
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccf_and_vars_chunklist(chunk_list::AbstractChunkList,
    plan_for_chunk::AbstractVector{PlanT} = BasicCCFPlan(); assume_sorted::Bool = false, use_pixel_vars::Bool = false  ) where {
                        PlanT<:AbstractCCFPlan }
    num_chunks = length(chunk_list)
    @assert length(plan_for_chunk) == num_chunks
    @assert assume_sorted || issorted( first(plan_for_chunk).line_list.λ )

    num_vels = maximum(map(chid->calc_length_ccf_v_grid(plan_for_chunk[chid]), 1:num_chunks ))
    chunk_ccfs_out = zeros( num_vels, num_chunks )
    chunk_ccf_vars_out = zeros( num_vels, num_chunks )

    calc_order_ccf_and_vars_chunklist!(chunk_ccfs_out, chunk_ccf_vars_out, chunk_list, plan_for_chunk, assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars)
    return (ccfs=chunk_ccfs_out, ccf_vars=chunk_ccf_vars_out)
end
