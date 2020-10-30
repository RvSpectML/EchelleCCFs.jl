"""
Functions for estimating a template CCF

Author: Eric Ford
Created: October 2020
"""

function calc_normalized_ccfs( ccfs::A2 ) where {T2<:Real, A2<:AbstractArray{T2,2} }
    return ccfs ./maximum(ccfs,dims=1)
end

function calc_normalized_ccfs( ccfs::A2, ccf_vars::A3 ) where { T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    max_ccfs = maximum(ccfs,dims=1)
    ccfs_norm = ccfs ./ max_ccfs
    ccf_vars_norm = ccf_vars ./ max_ccfs
    return (ccfs=ccfs_norm, ccf_vars=ccf_vars_norm)
end

function calc_ccf_template( ccfs::A2; assume_normalized::Bool = false ) where {T2<:Real, A2<:AbstractArray{T2,2} }
    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    ccf_template = mean(ccfs_norm,dims=2)
    est_ccf_vars_norm = var(ccfs_norm.-ccf_template,dims=2)
    ccf_template = vec(sum(ccfs_norm./est_ccf_vars_norm,dims=2) ./sum(1.0 ./est_ccf_vars_norm,dims=2))
    return ccf_template
end

function calc_ccf_template( ccfs::A2, ccf_vars::A3; assume_normalized::Bool = false  ) where { T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    (ccfs_norm, ccf_vars_norm) = assume_normalized ? (ccfs, ccf_vars) : calc_normalized_ccfs(ccfs, ccf_vars)
    ccf_template = vec(sum(ccfs./ccf_vars,dims=2) ./sum(1.0 ./ccf_vars,dims=2))
    return ccf_template
end
