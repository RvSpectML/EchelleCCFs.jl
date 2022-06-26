"""
Functions for estimating covariance matrices to use when fitting CCFs

Author: Eric Ford
Created: October 2020
"""

""" `calc_ccf_sample_covar( ccfs, ccf_vars)`
Estimate sample covariance for set of CCFs
"""
function calc_ccf_sample_covar( ccfs::A2, ccf_vars::A3; assume_normalized::Bool = false, verbose::Bool = false ) where { T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    num_vels, num_obs = size(ccfs)
    (ccfs_norm, ccf_vars_norm) = assume_normalized ? (ccfs, ccf_vars) : calc_normalized_ccfs(ccfs, ccf_vars)
    ccf_template = calc_ccf_template(ccfs_norm, ccf_vars_norm, assume_normalized=true )
    obs_weights = calc_ccf_weights_formal(ccfs, ccf_vars, ccf_template)
    if verbose   println("median(obs_weights) = ", median(obs_weights), "   extrema(obs_weights) = ", extrema(obs_weights))   end
    ccf_sample_covar = Array{eltype(ccfs),2}(undef, num_vels, num_vels)
    ccf_sample_covar = calc_ccf_sample_covar_helper!(ccf_sample_covar, ccfs_norm, obs_weights, template=ccf_template, assume_normalized=true )
    ccf_sample_covar = Symmetric(ccf_sample_covar)
    return ccf_sample_covar
end

#=
function calc_ccf_sample_covar( ccfs::A2; assume_normalized::Bool = false  ) where {T2<:Real, A2<:AbstractArray{T2,2} }
    @warn "You probably shouldn't be using this without also passing ccf_var"
    num_vels, num_obs = size(ccfs)
    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    ccf_template = calc_ccf_template(ccfs_norm, assume_normalized=true)
    obs_weights = calc_ccf_weights_empirical(ccfs_norm, ccf_template)
    println("median(obs_weights) = ", median(obs_weights), "   extrema(obs_weights) = ", extrema(obs_weights))
    ccf_sample_covar = Array{eltype(ccfs),2}(undef, num_vels, num_vels)
    calc_ccf_sample_covar_helper!(ccf_sample_covar, ccfs_norm, obs_weights, template=ccf_template, assume_normalized=true )
    return ccf_sample_covar
end
=#

function calc_ccf_sample_covar_helper!( covar_out::A2, ccfs::A2; template::A1, assume_normalized::Bool = false  ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    num_vels, num_obs = size(ccfs)
    @assert 5 <= num_vels <= 10000  # Arbitrary limits
    @assert 2 <= num_obs <= 10000   # Arbitrary upper limit
    @assert size(covar_out,1) == size(covar_out,2) == num_vels
    @assert length(template) == num_vels
    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    covar_out .= ((ccfs_norm .- template) * (ccfs_norm .- template)') ./ (num_vels-1)
    return covar_out
end

function calc_ccf_sample_covar_helper!( covar_out::A2, ccfs::A2, weights::A1; template::A1, assume_normalized::Bool = false  ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    num_vels, num_obs = size(ccfs)
    @assert 5 <= size(ccfs,1) <= 10000  # Arbitrary limits
    @assert 2 <= num_obs <= 10000       # Arbitrary upper limit
    @assert size(covar_out,1) == size(covar_out,2) == num_vels
    @assert length(weights) == num_obs
    @assert length(template) == num_vels
    @assert sum(weights) ≈ one(eltype(weights))

    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    # Weighted version of (ccfs_norm .- template) * (ccfs_norm .- template)'
    covar_out = zeros(num_vels, num_vels)
    for j in 1:num_vels
        for k in 1:j
            for i in 1:num_obs
                summand =  weights[i] * (ccfs[j,i]-template[j]) * (ccfs[k,i]-template[k])
                covar_out[k,j] += summand
                #=  Can skip since wrapping with Symmetric below to improve memory access pattern
                if j<k
                    covar_out[k,j] += summand
                end
                =#
            end # i
        end # j
    end # k
    sum_w = sum(weights)
    sum_w2 = sum(weights.^2)
    n_eff = sum_w^2/sum_w2
    normalization = 1.0 / (1.0 - sum_w2 )
    normalization *= n_eff/(n_eff-1)
    covar_out *= normalization
    return covar_out
end


# TODO: include weighting by template derivative, so SNR per velocity pixel is weighted by information content
function calc_ccf_weights_formal(ccfs::A2, ccf_vars::A3, template::A1; assume_normalized::Bool = false, verbose::Bool = false ) where {  T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    weights = vec(1.0 ./ sum(ccf_vars./ccfs,dims=1))
    weights ./= sum(weights)
end

# TODO: include weighting by template derivative, so SNR per velocity pixel is weighted by information content
# TODO: Could use difference relative to best-fitting shifted template CCF
# TODO: Could use better CCF model, e.g., reduced rank version
function calc_ccf_weights_empirical(ccfs::A2, template::A1; assume_normalized::Bool = false, verbose::Bool = false ) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    weights = vec(1.0 ./ sum((ccfs .- template).^2,dims=1))
    weights ./= sum(weights)
end
