"""
Functions for estimating covariance matrices to use when fitting CCFs

Author: Eric Ford
Created: October 2020
"""

default_kmax_calc_ccf_sample_covar_reduced_rank = 6  # TODO:  pick more thoughtfully than just using values from Scapels paper

""" `calc_ccf_sample_covar_reduced_rank( ccfs, ccf_vars ; kmax)`
Calculate a reduced rank approximation to the sample covariance matrix for the CCF
"""
function calc_ccf_sample_covar_reduced_rank( ccfs::A2, ccf_vars::A3 ; kmax::Integer = default_kmax_calc_ccf_sample_covar_reduced_rank, assume_normalized::Bool = false, verbose::Bool = false ) where {T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    num_obs = size(ccfs,2)
    (ccfs_norm, ccf_vars_norm) = assume_normalized ? (ccfs, ccf_vars) : calc_normalized_ccfs(ccfs, ccf_vars)
    ccf_template = calc_ccf_template(ccfs_norm, ccf_vars_norm, assume_normalized=true )
    return calc_ccf_sample_covar_reduced_rank_helper(ccfs_norm, template=ccf_template, kmax=kmax)
end

#=
function calc_ccf_sample_covar_reduced_rank( ccfs::A2; kmax::Integer = default_kmax_calc_ccf_sample_covar_reduced_rank, assume_normalized::Bool = false ) where {T2<:Real, A2<:AbstractArray{T2,2} }
    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    ccf_template = calc_ccf_template(ccfs_norm, assume_normalized=true )
    #ccf_sample_covar = (ccfs_norm .- template) * (ccfs_norm .- template)'
    return calc_ccf_sample_covar_reduced_rank_helper(ccfs_norm, template=ccf_template, kmax=kmax)
end
=#

function calc_ccf_sample_covar_reduced_rank_helper( ccfs_norm::A2; template::A1, kmax::Integer = default_kmax_calc_ccf_sample_covar_reduced_rank, verbose::Bool = false ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    num_obs = size(ccfs_norm,2)
    @assert 2 <= num_obs <= 10000  # Arbitrary upper limit
    sample_covar_svd = svd((ccfs_norm .- template)')
    R_kmax = sample_covar_svd.U[1:kmax,:]'*(diagm(sample_covar_svd.S[1:kmax])*sample_covar_svd.Vt[1:kmax,:])
    #ccf_sample_covar = (max.(R_kmax'*R_kmax,0)) / num_obs
    ccf_sample_covar = (R_kmax'*R_kmax)  / (num_obs-1)
    return (covar=ccf_sample_covar, R_kmax=R_kmax)
end

""" `est_covar_for_obs(covar_reduced_rank, diag_model, diag_norm)`
Returns estimate of the CCF covariance based using
`covar_reduced_rank + diag_norm * diag_model`
`diag_model` is typically not diagonal, but rather terms near the diagonal
"""
function est_covar_for_obs(covar_reduced_rank::A2, diag_model::A3, diag_norm::Real) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    @assert size(covar_reduced_rank) == size(diag_model)
    return covar_reduced_rank + diag_norm * diag_model
end

function est_covar_for_obs(covar_reduced_rank::A2, shape_model::A3, row_norms::A1, obs_idx::Integer) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    @assert 1 <= obs_idx <= length(row_norms)
    est_covar_for_obs(covar_reduced_rank, shape_model, row_norms[obs_idx])
end


@. triangle_helper(x, p) = p[1]*max(1-abs(x)/p[2], 0)

function triangular_shape_near_diagonal(n::Integer, width::Real; norm::Real = 1)
  shape_model = collect( triangle_helper(i-j,[norm,width]) for i in 1:n, j in 1:n )
end


#=
using LsqFit

# WARNING:  This function still needs major work
function calc_ccf_sample_covar_and_near_diag_model( ccfs::A2, ccf_vars::A3 ; kmax::Integer = default_kmax_calc_ccf_sample_covar_reduced_rank,
                                            num_diags_to_calc::Integer = floor(Int,size(ccfs,1)//2), verbose::Bool = false ) where {T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    num_vels, num_obs = size(ccfs)
    @assert size(ccfs) == size(ccf_vars)
    @assert 5 <= num_diags_to_calc <= num_vels
    @assert kmax <= num_obs
    (ccfs_norm, ccf_vars_norm) = #= assume_normalized ? (ccfs, ccf_vars) : =# calc_normalized_ccfs(ccfs, ccf_vars)
    ccf_template = calc_ccf_template(ccfs_norm, ccf_vars_norm, assume_normalized=true )
    obs_weights = calc_ccf_weights_formal(ccfs, ccf_vars, ccf_template)
    ccf_sample_covar = Array{eltype(ccfs),2}(undef, num_vels, num_vels)
    calc_ccf_sample_covar_helper!(ccf_sample_covar, ccfs_norm, obs_weights, template=ccf_template, assume_normalized=true )

    (ccf_sample_covar_reduced_rank, R_kmax) = calc_ccf_sample_covar_reduced_rank_helper( ccfs_norm, template = ccf_template, kmax=kmax) #, assume_normalized=true )
    residR = ccf_sample_covar - ccf_sample_covar_reduced_rank
    return residR
    if verbose
        flush(stdout)
        println("# extrema(ccf_sample_covar) = ",extrema(ccf_sample_covar))
        println("# extrema(ccf_sample_covar_reduced_rank) = ",extrema(ccf_sample_covar_reduced_rank))
        println("# extrema(residR) = ",extrema(residR))
    end
    #=
    mean_zero_lag = compute_mean_off_diag(residR,0)
    mean_offset_diag = map(i->compute_mean_off_diag(residR,i)./mean_zero_lag, -num_diags_to_calc:num_diags_to_calc)
    =#
    # TODO: fit triangle away from CCF peak
    residR_away_from_peak_1 = view(residR,1:ceil(Int,num_vels//4),1:ceil(Int,num_vels//4))
    residR_away_from_peak_2 = view(residR,floor(Int,num_vels*3//4):num_vels,floor(Int,num_vels*3//4):num_vels)
    mean_zero_lag = (compute_mean_off_diag(residR_away_from_peak_1,0)+compute_mean_off_diag(residR_away_from_peak_2,0))/2
    println("mean_zero_lag = ",mean_zero_lag)
    num_diags_to_calc_tmp_1 = min(num_diags_to_calc,size(residR_away_from_peak_1,1)-1)
    num_diags_to_calc_tmp_2 = min(num_diags_to_calc,size(residR_away_from_peak_2,1)-1)
    num_diags_to_calc_tmp = min(num_diags_to_calc_tmp_1,num_diags_to_calc_tmp_2)
    println("num_diags_to_calc_tmp = ",num_diags_to_calc_tmp)
    mean_offset_diag   = map(i->compute_mean_off_diag(residR_away_from_peak_1,i)./mean_zero_lag, -num_diags_to_calc_tmp:num_diags_to_calc_tmp)
    println("mean_offset_diag = ",extrema(mean_offset_diag) )
    mean_offset_diag .+= map(i->compute_mean_off_diag(residR_away_from_peak_2,i)./mean_zero_lag, -num_diags_to_calc_tmp:num_diags_to_calc_tmp)
    mean_offset_diag ./= 2
    println("mean_offset_diag = ",extrema(mean_offset_diag) )
    return mean_offset_diag
    # Find idx around CCF peak
    #mid_offset = findmin(ccf_template)[2]
    half_width = est_full_width(1:num_vels,ccf_template)
    width_guess = 2*half_width
    result = curve_fit(triangle_helper, -num_diags_to_calc_tmp:num_diags_to_calc_tmp, mean_offset_diag, [1.0, width_guess] )
    #=
    max_offset = min(num_vels,mid_offset+num_diags_to_calc)
    min_offset = max(1,mid_offset-num_diags_to_calc)
    epsilon = 0.1 * mean_zero_lag
    idx_hi = findfirst(x->x<=epsilon,view(mean_offset_diag,mid_offset:max_offset))
    idx_lo = findfirst(x->x<=epsilon,view(mean_offset_diag,mid_offset:-1:min_offset))
    if idx_hi == nothing || idx_lo == nothing
        #flush(stdout)
        #println("mean_offset_diag = ",view(mean_offset_diag,mid_offset:max_offset))
        idx_hi = min(num_vels,mid_offset + 2*half_width)
        idx_lo = max(1,mid_offset - 2*half_width)
    else
        idx_hi = min(num_vels, idx_hi+ (mid_offset-1) )
        idx_lo = max(1, -idx_lo + (mid_offset+1) )
    end
    width_guess = idx_hi-idx_lo
    if verbose   println("# idx_lo = ", idx_lo, " idx_hi = ", idx_hi, " width_guess = ", width_guess )   end
    result = curve_fit(triangle_helper, view(-num_diags_to_calc:num_diags_to_calc,idx_lo:idx_hi), view(mean_offset_diag,idx_lo:idx_hi), [1.0, width_guess] )
    =#
    A = coef(result)[1]
    δv = coef(result)[2]
    if true || verbose   println("# A = ", A, " δv = ", δv, "  mean_lag_zero = ", mean_zero_lag, " converged = ", result.converged )   end
    #shape_model = ( triangle_helper(i-j,[A,δv]) for i in 1:num_vels, j in 1:num_vels )
    shape_model = triangular_shape_near_diagonal(num_vels, δv; norm=A)
    row_norms = vec(var( (ccfs_norm .- ccf_template)' .-  R_kmax, dims=2))
    #return (covar_reduced_rank=ccf_sample_covar_reduced_rank / num_obs, shape_model=shape_model, row_norms=row_norms)
    return (covar_reduced_rank=ccf_sample_covar_reduced_rank , shape_model=shape_model, row_norms=row_norms)
end
=#

""" ` compute_mean_off_diag(matrix, offset)`
Compute mean value of entries with given offset from the diagonal.
Note: offset may be positive or negative
"""
function compute_mean_off_diag(m::AbstractArray{T,2}, offset::Integer) where {T<:Real}
  n = size(m,1)
  @assert size(m,1) == size(m,2)
  @assert -n <= offset <= n
  sum = 0
  count = 0
  for i in 1:n
      if 1<= i+offset <= n
        sum += m[i,i+offset]
        count += 1
      end
  end
  return sum/count
end
