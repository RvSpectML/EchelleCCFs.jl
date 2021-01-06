function assign_lines_to_orders(line_list::DataFrame, order_info::DataFrame; v_center::Real=0 )
  ll1 = copy(line_list)
  ll2 = copy(line_list)
  ll1[!,:order] .= 0
  ll2[!,:order] .= 0
  for (i,line) in enumerate(eachrow(line_list))
    order_indices = findall( map(order-> order.λ_min <= line.lambda_lo <= order.λ_max && order.λ_min <= line.lambda_hi <= order.λ_max, eachrow(order_info) ) )
    if length(order_indices) == 1
      ll1[i,:order] = order_info.order[order_indices[1]]
    elseif length(order_indices) == 2
      ll1[i,:order] = order_info.order[order_indices[1]]
      ll2[i,:order] = order_info.order[order_indices[2]]
    elseif length(order_indices) > 2
      @warn "One line appears in >2 orders.  Not implemented yet."
      println("line = ", line, "  order_indices = ", order_indices, "  line.order = ",order_info.order[order_indices])
      #=  df_tmp = DataFrame()
      for key in keys(eachcol(line_list))
        df_tmp[!,key] = fill(line[key], length(order_indices)-1 )
      end
      df_tmp[:order] .= order_indices[2:end]
      append!(line_list,df_tmp)
      =#
      #break
    end
  end
  append!(ll1,ll2)
  ll1 |> @filter( _.order >= 1 ) |> @orderby(_.order) |> @thenby(_.lambda_lo) |> DataFrame
  #ll1 |> @filter( _.order >= 1 ) |> @orderby(_.lambda_lo) |> @thenby(_.order) |> DataFrame
end


function calc_snr_weights_for_lines!(line_list::DataFrame, spectra::AAS ; v_center::Real=0 ) where { AS<:AbstractSpectra2D, AAS<:AbstractVector{AS} }
  num_obs = length(spectra)
  @assert num_obs >= 1
  num_cols = size(first(spectra).λ,1)
  num_lines = size(line_list,1)
  snr_list = zeros(num_obs, num_lines)
  depth = zeros(num_obs, num_lines)
  exp_sigma_rv = zeros(num_obs, num_lines)
  line_list[!,:weight_snr] .= 0.0
  for (i,line) in enumerate(eachrow(line_list))
    order = line.order
    for obs_idx in 1:num_obs
      boost_factor = calc_doppler_factor(v_center)
      # TODO: Do something with boost_factor for each observation
      idx_min = searchsortedfirst(view(spectra[obs_idx].λ,1:num_cols,order),line.lambda_lo)
      idx_max = searchsortedlast(view(spectra[obs_idx].λ,1:num_cols,order),line.lambda_hi)
      @assert 1 <= idx_min <= num_cols
      @assert 1 <= idx_max <= num_cols
      #println("# Line ", i, " at ", line.lambda, " obs ", obs_idx, " idx = ", idx_min, " - ", idx_max)
      snr_list[obs_idx,i] = RvSpectMLBase.calc_snr( spectra[obs_idx], idx_min:idx_max, order )
      #(depth[obs_idx,i], exp_sigma_rv[obs_idx,i]) = calc_depth_and_expected_rv_precission( spectra[obs_idx], idx_min:idx_max, order )
    end
  end
  snr_weights = vec(median(snr_list,dims=1)).^2
  sum_snr_weights = sum(snr_weights)
  sum_weight_inputs = sum(line_list[!,:weight])
  snr_weights .*= sum_weight_inputs/sum_snr_weights
  line_list[!,:weight_snr] .= snr_weights
  #=
  line_list[!,:exp_sigma_rv] .= vec(median(exp_sigma_rv,dims=1))
  exp_σ_rv_weights = 1.0 ./ line_list[!,:exp_sigma_rv].^2
  sum_σ_rv_weights = sum(exp_σ_rv_weights)
  exp_σ_rv_weights ./= sum_σ_rv_weights

  line_list[!,:weight_theory] .= exp_σ_rv_weights #.* snr_weights
  #line_list[!,:weight_theory] ./= sum(line_list[!,:weight_theory])
  =#

  line_list[!,:weight_input] = deepcopy(line_list[!,:weight])
  #line_list[!,:weight] = line_list[!,:weight_theory]
  line_list[!,:weight] .*= snr_weights
  #line_list[!,:depth] .= vec(median(depth,dims=1))
  return line_list
end
