"""
    Code for convenience functions for calculating CCFs for each chunk of an AbstractChunkListTimeseries (see RvSpectMLBase)
Author: Eric Ford
Created: August 2020
"""


"""  `calc_order_ccf_chunklist_timeseries( chunklist_timeseries, ccf_plan )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) of each spectrum in a timeseries.
    CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
# Return:
A 3-d array containing the CCF at each (velocity, chunk, spectrum)
Note that the ccf_plan provided is used as a template for creating a custom ccf_plan for each chunk that
    only includes lines that reliably appear in that order for all spectra in the chunklist_timeseries.
"""
function calc_order_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                plan::PlanT = BasicCCFPlan(); verbose::Bool = false ,
                Δfwhm::AbstractVector{T} = zeros(0)
                ) where { PlanT<:AbstractCCFPlan, T<:Real  }

  #@assert issorted( plan.line_list.λ )
  #nvs = length(calc_ccf_v_grid(plan))
  nvs = calc_length_ccf_v_grid(plan)
  norders = num_chunks(clt) # length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nvs, norders, nobs)
  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,norders)
  for chid in 1:norders
      order = clt[1].order[chid]
      if typeof(plan.line_list) <: BasicLineList2D
          start_order_idx = searchsortedfirst(plan.line_list.order,order)
          stop_order_idx = searchsortedlast(view(plan.line_list.order,start_order_idx:num_lines),order) + start_order_idx-1
      else
          start_order_idx = 1
          stop_order_idx = num_lines
      end
      if (1<=start_order_idx<=stop_order_idx<=num_lines)
          # find the maximum lower wavelength for the chunk, and the minumum upper wavelength, over all observations
          λmin = maximum(map(obsid->first(clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
          λmax = minimum(map(obsid->last( clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
          # extend λmin/λmax by the velocity range over which we don't want the mask to change
          λmin  = λmin / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(-plan.v_range_no_mask_change))
          λmax  = λmax / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(plan.v_range_no_mask_change))
          # extend λmin/λmax by the mask width
          upper_edge_of_mask_for_line_at_λmin = λ_max(plan.mask_shape,λmin)
          lower_edge_of_mask_for_line_at_λmax = λ_min(plan.mask_shape,λmax)
          # find the first and last mask entries to use in each chunk
          start_line_idx = searchsortedfirst(view(plan.line_list.λ,start_order_idx:stop_order_idx),upper_edge_of_mask_for_line_at_λmin) + start_order_idx-1
          #stop_line_idx  = num_lines+1 - searchsortedfirst(view(plan.line_list.λ,num_lines:-1:1),lower_edge_of_mask_for_line_at_λmax,rev=true)
          stop_line_idx  = start_order_idx-1 + searchsortedlast(view(plan.line_list.λ,start_order_idx:stop_order_idx),lower_edge_of_mask_for_line_at_λmax)
          if verbose
              flush(stdout)
              println("extrema(plan.line_list.λ) = ",extrema(plan.line_list.λ) )
              if (1 <= start_line_idx <= length(plan.line_list.λ)) && (1 <= stop_line_idx <= length(plan.line_list.λ))
                  println("start_order_idx = ", start_order_idx, "  stop_order_idx = ", stop_order_idx)
                  println("start_line_idx = ", start_line_idx, " λ= ", plan.line_list.λ[start_line_idx]) #, " order= ", plan.line_list.order[start_line_idx])
                  println("stop_line_idx = ", stop_line_idx, " λ= ", plan.line_list.λ[stop_line_idx]) #, " order= ", plan.line_list.order[stop_line_idx])
              end
              println("upper_edge_of_mask_for_line_at_λmin = ", upper_edge_of_mask_for_line_at_λmin, " Lower_edge_of_mask_for_line_at_λmax = ", lower_edge_of_mask_for_line_at_λmax)
          end

          if (1 <= start_line_idx <= length(plan.line_list.λ)) && (1 <= stop_line_idx <= length(plan.line_list.λ))
              if verbose && typeof(plan.line_list) <: BasicLineList2D
                  println("start_line_idx = ", start_line_idx, " λ= ", plan.line_list.λ[start_line_idx], " order= ", plan.line_list.order[start_line_idx])
                  println("stop_line_idx = ", stop_line_idx, " λ= ", plan.line_list.λ[stop_line_idx], " order= ", plan.line_list.order[stop_line_idx])
                  println(" Using ", length(start_line_idx:stop_line_idx), " lines for chunk ", chid, " order ", order)
              end
              if ! issorted( view(plan.line_list.λ, start_line_idx:stop_line_idx ) )
                  println("# View is not sorted.  Why?")
                  for i in start_line_idx:stop_line_idx
                      println("i= ", i, " λ= ", plan.line_list.λ[i], " order= ", plan.line_list.order[i])
                  end
              end
              @assert issorted( view(plan.line_list.λ, start_line_idx:stop_line_idx ) )
              # line_list_for_chunk = BasicLineList2D(view(plan.line_list.λ,start_line_idx:stop_line_idx), view(plan.line_list.weight,start_line_idx:stop_line_idx), fill(order,length(start_line_idx:stop_line_idx)) )
              line_list_for_chunk = BasicLineList(view(plan.line_list.λ,start_line_idx:stop_line_idx), view(plan.line_list.weight,start_line_idx:stop_line_idx)  )
          else  # No lines in this chunk!
              line_list_for_chunk = EmptyBasicLineList()
          end
      else
        line_list_for_chunk = EmptyBasicLineList()
      end

      #create a plan for this chunk that only includes the mask entries we want for this chunk
      plan_for_chunk[chid] = BasicCCFPlan( line_list=line_list_for_chunk, midpoint=plan.v_center, step=plan.v_step, max=plan.v_max, mask_shape=plan.mask_shape )
  end
  #@threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk),hcat, 1:length(clt) )
  #Threads.@threads for i in 1:nobs
    #  order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan)
  #end
  Threads.@threads for i in 1:nobs
      if verbose 
         @info "Computing CCF: " i nobs (tid=Threads.threadid())
      end 
      this_Δfwhm = length(Δfwhm) == nobs ? Δfwhm[i] : 0.0
          max_nan_frac_to_use = 0.75
          nan_frac_lambda = sum(isnan.(clt.chunk_list[i].data[j].λ)) / length(clt.chunk_list[i].data[j].λ)
          nan_frac_flux =   sum(isnan.(clt.chunk_list[i].data[j].flux)) / length(clt.chunk_list[i].data[j].flux)
          nan_frac_var =    sum(isnan.(clt.chunk_list[i].data[j].var)) / length(clt.chunk_list[i].data[j].var)
          if (nan_frac_lambda > max_nan_frac_to_use) || (nan_frac_flux > max_nan_frac_to_use) || (nan_frac_var > max_nan_frac_to_use)  
             continue
          end 
      order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan_for_chunk, Δfwhm=this_Δfwhm, assume_sorted=true )
  end

  return order_ccfs
end

"""  `calc_order_ccf_and_var_chunklist_timeseries( chunklist_timeseries, ccf_plan )`
Convenience function to compute separate CCFs for each chunk (potentially an order or view around one or two lines) of each spectrum in a timeseries.
    CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
# Return:
A 3-d array containing the CCF at each (velocity, chunk, spectrum)
Note that the ccf_plan provided is used as a template for creating a custom ccf_plan for each chunk that
    only includes lines that reliably appear in that order for all spectra in the chunklist_timeseries.
"""
function calc_order_ccf_and_var_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                plan::PlanT = BasicCCFPlan(); verbose::Bool = false,
                ccf_var_scale::Real = 1.0, Δfwhm::AbstractVector{T} = zeros(0)
                 ) where {
                    PlanT<:AbstractCCFPlan, T<:Real }
  #@assert issorted( plan.line_list.λ )
  #nvs = length(calc_ccf_v_grid(plan))
  nvs = calc_length_ccf_v_grid(plan)
  norders = num_chunks(clt) # length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nvs, norders, nobs)
  order_ccf_vars = zeros(nvs, norders, nobs)
  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,norders)
  for chid in 1:norders
      order = clt[1].order[chid]
      if typeof(plan.line_list) <: BasicLineList2D
          start_order_idx = searchsortedfirst(plan.line_list.order,order)
          stop_order_idx = searchsortedlast(view(plan.line_list.order,start_order_idx:num_lines),order) + start_order_idx-1
      else
          start_order_idx = 1
          stop_order_idx = num_lines
      end
      if (1<=start_order_idx<=stop_order_idx<=num_lines)
          # find the maximum lower wavelength for the chunk, and the minumum upper wavelength, over all observations
          λmin = maximum(map(obsid->first(clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
          λmax = minimum(map(obsid->last( clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
          # extend λmin/λmax by the velocity range over which we don't want the mask to change
          λmin  = λmin / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(-plan.v_range_no_mask_change))
          λmax  = λmax / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(plan.v_range_no_mask_change))
          # extend λmin/λmax by the mask width
          upper_edge_of_mask_for_line_at_λmin = λ_max(plan.mask_shape,λmin)
          lower_edge_of_mask_for_line_at_λmax = λ_min(plan.mask_shape,λmax)
          # find the first and last mask entries to use in each chunk
          start_line_idx = searchsortedfirst(view(plan.line_list.λ,start_order_idx:stop_order_idx),upper_edge_of_mask_for_line_at_λmin) + start_order_idx-1
          stop_line_idx  = start_order_idx-1 + searchsortedlast(view(plan.line_list.λ,start_order_idx:stop_order_idx),lower_edge_of_mask_for_line_at_λmax)
          if verbose
              flush(stdout)
              println("extrema(plan.line_list.λ) = ",extrema(plan.line_list.λ) )
              if (1 <= start_line_idx <= length(plan.line_list.λ)) && (1 <= stop_line_idx <= length(plan.line_list.λ))
                  println("start_order_idx = ", start_order_idx, "  stop_order_idx = ", stop_order_idx)
                  println("start_line_idx = ", start_line_idx, " λ= ", plan.line_list.λ[start_line_idx]) #, " order= ", plan.line_list.order[start_line_idx])
                  println("stop_line_idx = ", stop_line_idx, " λ= ", plan.line_list.λ[stop_line_idx]) #, " order= ", plan.line_list.order[stop_line_idx])
              end
              println("upper_edge_of_mask_for_line_at_λmin = ", upper_edge_of_mask_for_line_at_λmin, " Lower_edge_of_mask_for_line_at_λmax = ", lower_edge_of_mask_for_line_at_λmax)
          end

          if (1 <= start_line_idx <= length(plan.line_list.λ)) && (1 <= stop_line_idx <= length(plan.line_list.λ))
              if verbose && typeof(plan.line_list) <: BasicLineList2D
                  println("start_line_idx = ", start_line_idx, " λ= ", plan.line_list.λ[start_line_idx], " order= ", plan.line_list.order[start_line_idx])
                  println("stop_line_idx = ", stop_line_idx, " λ= ", plan.line_list.λ[stop_line_idx], " order= ", plan.line_list.order[stop_line_idx])
                  println(" Using ", length(start_line_idx:stop_line_idx), " lines for chunk ", chid, " order ", order)
              end
              if ! issorted( view(plan.line_list.λ, start_line_idx:stop_line_idx ) )
                  println("# View is not sorted.  Why?")
                  for i in start_line_idx:stop_line_idx
                      println("i= ", i, " λ= ", plan.line_list.λ[i], " order= ", plan.line_list.order[i])
                  end
              end
              @assert issorted( view(plan.line_list.λ, start_line_idx:stop_line_idx ) )
              # line_list_for_chunk = BasicLineList(view(plan.line_list.λ,start_line_idx:stop_line_idx), view(plan.line_list.weight,start_line_idx:stop_line_idx) )
              line_list_for_chunk = BasicLineList(view(plan.line_list.λ,start_line_idx:stop_line_idx), view(plan.line_list.weight,start_line_idx:stop_line_idx)  )
          else  # No lines in this chunk!
              line_list_for_chunk = EmptyBasicLineList()
          end
      else
        line_list_for_chunk = EmptyBasicLineList()
      end
      #create a plan for this chunk that only includes the mask entries we want for this chunk
      plan_for_chunk[chid] = BasicCCFPlan( line_list=line_list_for_chunk, midpoint=plan.v_center, step=plan.v_step, max=plan.v_max, mask_shape=plan.mask_shape )
  end
  #@threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk),hcat, 1:length(clt) )
  #Threads.@threads for i in 1:nobs
    #  order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan)
  #end
  Threads.@threads for i in 1:nobs
      if verbose 
         @info "Computing CCF: " i nobs (tid=Threads.threadid())
      end
      this_Δfwhm = length(Δfwhm) == nobs ? Δfwhm[i] : 0.0
      max_nan_frac_to_use = 0.75
      for j in 1:norders
          nan_frac_lambda = sum(isnan.(clt.chunk_list[i].data[j].λ)) / length(clt.chunk_list[i].data[j].λ)
          nan_frac_flux =   sum(isnan.(clt.chunk_list[i].data[j].flux)) / length(clt.chunk_list[i].data[j].flux)
          nan_frac_var =    sum(isnan.(clt.chunk_list[i].data[j].var)) / length(clt.chunk_list[i].data[j].var)
          if (nan_frac_lambda > max_nan_frac_to_use) || (nan_frac_flux > max_nan_frac_to_use) || (nan_frac_var > max_nan_frac_to_use)  
             continue
          end 
          ( ccf_tmp, ccf_var_tmp ) = calc_ccf_and_var_chunk(clt.chunk_list[i][j], plan_for_chunk[j],
                                            ccf_var_scale=ccf_var_scale, Δfwhm=this_Δfwhm, assume_sorted=true )
          order_ccfs[:,j,i] .= ccf_tmp
          order_ccf_vars[:,j,i] .= ccf_var_tmp
      end
  end

  return (ccf=order_ccfs, ccf_var=order_ccf_vars)
end
