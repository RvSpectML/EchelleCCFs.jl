"""
    Code for convenience functions for calculating summed CCFs for each spectra in an AbstractChunkListTimeseries (see RvSpectMLBase)
Author: Eric Ford
Created: August 2020
"""

"""  `calc_ccf_chunklist_timeseries( chunklist_timeseries, line_list )`
Convenience function to compute CCF for a timeseries of spectra, each with a chunklist.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
- verbose (false)
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the ccf_plan.
Note that the ccf_plan provided is used as a template for creating a custom ccf_plan for each chunk that
    only includes lines that reliably appear in that order for all spectra in the chunklist_timeseries.
"""
function calc_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                                plan::PlanT = BasicCCFPlan(); verbose::Bool = false,
                                calc_ccf_var::Bool = false, use_pixel_vars::Bool = false  ) where { PlanT<:AbstractCCFPlan }

  @assert length(clt) >= 1
  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,num_chunks(clt))
  if use_pixel_vars
      var_list = Vector{Vector{Float64}}(undef,num_chunks(clt))
  end
  for chid in 1:num_chunks(clt)
      order = clt[1].order[chid]
      if typeof(plan.line_list) <: BasicLineList2D
          start_order_idx = searchsortedfirst(plan.line_list.order,order)
          stop_order_idx = searchsortedlast(view(plan.line_list.order,start_order_idx:num_lines),order) + start_order_idx-1
      else
          start_order_idx = 1
          stop_order_idx = num_lines
      end

      if verbose
          println("# Chunk ", chid, " order= ", order)
          println("# start_order_idx = ", start_order_idx, "    stop_order_idx = ", stop_order_idx)
              if (1<=start_order_idx<=stop_order_idx<=num_lines) && typeof(plan.line_list) <: BasicLineList2D
                  println("   order = ", plan.line_list.order[start_order_idx], " - ", plan.line_list.order[stop_order_idx])
              end
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
      if use_pixel_vars
          @warn "calc_ccf_chunklist_timeseries using pixel_vars is not yet implemented/tested."
          # Warning: Orders have consistent length over time, but chunks may not.  So currently only works with full orders
          var_list[chid] = mapreduce(obsid->clt.chunk_list[obsid].data[chid].var,+, 1:length(clt) )
          var_list[chid] ./= length(clt)
          #var_list[chid] = vec(median(mapreduce(obsid->clt.chunk_list[obsid].data[chid].var,hcat, 1:length(clt) ),dims=2))
            #return var_list
      end
  end
  #=
  if use_pixel_vars
      return @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], var_list, plan_for_chunk,assume_sorted=true, use_pixel_vars=true),hcat, 1:length(clt) )
  else
      return @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk,assume_sorted=true, use_pixel_vars=false),hcat, 1:length(clt) )
  end
  =#
  if calc_ccf_var
       #return @threaded mapreduce(obsid->calc_ccf_and_var_chunklist(clt.chunk_list[obsid], plan_for_chunk,assume_sorted=true, use_pixel_vars=true),hcat, 1:length(clt) )
       list_of_ccf_and_var = @threaded map(obsid->calc_ccf_and_var_chunklist(clt.chunk_list[obsid], plan_for_chunk,assume_sorted=true #=, use_pixel_vars=true =#), 1:length(clt) )
       nvs = length(first(list_of_ccf_and_var).ccf)
       ccfs_out = zeros(nvs,length(clt))
       ccf_vars_out = zeros(nvs,length(clt))
       for obsid in 1:length(clt)
           ccfs_out[:,obsid] .= list_of_ccf_and_var[obsid].ccf
           ccf_vars_out[:,obsid] .= list_of_ccf_and_var[obsid].ccf_var
       end
       return (ccfs=ccfs_out, ccf_vars=ccf_vars_out)
  else
       return @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk,assume_sorted=true #=, use_pixel_vars=false=# ),hcat, 1:length(clt) )
  end
end
