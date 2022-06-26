"""
Functions for estimating a template CCF

Author: Eric Ford
Created: October 2020
"""

""" calc_normalized_ccfs( ccfs )
Normalizes each spectrum by fitting a line to the region outside the center
"""
function fit_ccf_normalizations( v_grid::A1, ccfs::A2; v_mid::Real, dv_min::Real = RvSpectMLBase.max_bc, dv_max::Real = Inf ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    @assert zero(dv_min) < dv_min < Inf
    @assert dv_min < dv_max
    idx = findall( x->dv_min < abs(x-v_mid) < dv_max , v_grid )
    @assert length(idx) >= 2
    map(col-> mean(ccfs[idx,col] ), 1:size(ccfs,2) )
end

function calc_normalized_ccfs( v_grid::A1, ccfs::A2; v_mid::Real, dv_min::Real = RvSpectMLBase.max_bc, dv_max::Real = Inf ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    norms = fit_ccf_normalizations(v_grid,ccfs,v_mid=v_mid,dv_min=dv_min, dv_max=dv_max)
    ccfs_out = copy(ccfs)
    for i in 1:size(ccfs,2)
        ccfs_out[:,i] ./= norms[i]
    end
    return ccfs_out
end

function calc_normalized_ccfs( v_grid::A1, ccfs::A2, ccf_vars::A3; v_mid::Real, dv_min::Real = RvSpectMLBase.max_bc, dv_max::Real = Inf ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    norms = fit_ccf_normalizations(v_grid,ccfs,v_mid=v_mid,dv_min=dv_min, dv_max=dv_max)
    ccfs_out = copy(ccfs)
    ccf_vars_out = copy(ccf_vars)
    for i in 1:size(ccfs,2)
        ccfs_out[:,i] ./= norms[i]
        ccf_vars_out[:,i] ./= norms[i]^2
    end
    return (ccfs=ccfs_out, ccf_vars=ccf_vars_out)
end

""" calc_normalized_ccfs( ccfs )
Normalizes each spectrum by its maximum value.
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

""" calc_ccf_template( ccfs, [ccf_vars] ; assume_normalized = false )
Calculates ccf template
Warning: uses maximum CCF for normalization, unless you normalize manually.
"""
function calc_ccf_template( ccfs::A2; assume_normalized::Bool = false ) where {T2<:Real, A2<:AbstractArray{T2,2} }
    ccfs_norm = assume_normalized ? ccfs : calc_normalized_ccfs(ccfs)
    ccf_template = mean(ccfs_norm,dims=2)
    est_ccf_vars_norm = var(ccfs_norm.-ccf_template,dims=1)
    ccf_template = vec(sum(ccfs_norm./est_ccf_vars_norm,dims=2) ./sum(1.0 ./est_ccf_vars_norm,dims=2))
    return ccf_template
end

function calc_ccf_template( ccfs::A2, ccf_vars::A3; assume_normalized::Bool = false  ) where { T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    (ccfs_norm, ccf_vars_norm) = assume_normalized ? (ccfs, ccf_vars) : calc_normalized_ccfs(ccfs, ccf_vars)
    ccf_template = vec(sum(ccfs_norm./ccf_vars_norm,dims=2) ./sum(1.0 ./ccf_vars_norm,dims=2))
    return ccf_template
end
