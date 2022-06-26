"""
    Code to compute CCFs
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored and optimized by Eric Ford
"""

"""
    `ccf_1D!(ccf_out, λs, fluxes; ccf_plan )`

Compute the cross correlation function of a spectrum with a mask.
    Generalized version that should work with different mask shapes.
# Inputs:
- `ccf_out`: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- `λs`: 1-d array of wavelengths
- `fluxes`:  1-d array of fluxes
# Optional Arguments:
- `plan`:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D!(ccf_out::A1, λ::A2, flux::A3,
                plan::PlanT = BasicCCFPlan(); projection_workspace::AbstractArray{T4,2} = zeros(length(λ),1),
                assume_sorted::Bool = false ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1}, T4<:Real,
                PlanT<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    @assert assume_sorted || issorted( plan.line_list.λ )
    v_grid = calc_ccf_v_grid(plan)
    @assert size(ccf_out,1) == length(v_grid)

    if length(plan.line_list) < 1  # If no lines in chunk
        ccf_out .= zero(eltype(A1))
        return ccf_out
    end

    @assert size(projection_workspace,1) == length(λ)

    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:size(ccf_out,1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection_workspace, λ, plan, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = plan.allow_nans ? NaNMath.sum(flux .* projection_workspace) : sum(flux .* projection_workspace)
    end
    return ccf_out
end

"""
    `ccf_1D!(ccf_out, ccf_var_out, λs, fluxes, var; ccf_plan )`

Compute the cross correlation function of a spectrum with a mask.
    Generalized version that should work with different mask shapes.
# Inputs:
- `ccf_out`: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- `ccf_var_out`:  1-d array of size length(calc_ccf_v_grid(plan)) to store output
- `λs`: 1-d array of wavelengths
- `fluxes`:  1-d array of fluxes
- `var`:  1-d array of flux variances
# Optional Arguments:
- `plan`:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D!(ccf_out::A1, ccf_var_out::A1, λ::A2, flux::A3, var::A4,
                plan::PlanT = BasicCCFPlan(); projection_workspace::AbstractArray{T5,2} = zeros(length(λ),1),
                ccf_var_scale::Real = 1.0, assume_sorted::Bool = false ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1}, T4<:Real, A4<:AbstractArray{T4,1}, T5<:Real,
                PlanT<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    @assert assume_sorted || issorted( plan.line_list.λ )
    v_grid = calc_ccf_v_grid(plan)
    @assert size(ccf_out,1) == size(ccf_var_out,1)
    @assert size(ccf_out,1) == length(v_grid)


    if length(plan.line_list) < 1  # If no lines in chunk
        ccf_out .= zero(eltype(A1))
        ccf_var_out .= zero(eltype(A1))
        return (ccf=ccf_out, ccf_var=ccf_var_out)
    end

    @assert size(projection_workspace,1) == length(λ)
    @assert size(projection_workspace,2) == 1   # Keep until implement CCF w/ multiple masks simultaneously

    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:size(ccf_out,1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection_workspace, λ, plan, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = plan.allow_nans ? NaNMath.sum(flux .* projection_workspace) : sum(flux .* projection_workspace)
        ccf_var_out[i] = plan.allow_nans ? NaNMath.sum(var .* projection_workspace.^2) : sum(var .* projection_workspace.^2)
        ccf_var_out[i] *= ccf_var_scale
    end
    return (ccf=ccf_out, ccf_var=ccf_var_out)
end


"""
    `ccf_1D!(ccf_out, ccf_covar_out, λs, fluxes, var; ccf_plan )`

Compute the cross correlation function of a spectrum with a mask.
    Generalized version that should work with different mask shapes.
    WIP:  Generalized version to compute covariance matrix for ccf and optionally to account for LSF
# Inputs:
- `ccf_out`: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- `ccf_covar_out`:  2-d array of size length(calc_ccf_v_grid(plan))^2 to store output
- `λs`: 1-d array of wavelengths
- `fluxes`:  1-d array of fluxes
- `var`:  1-d array of flux variances
# Optional Arguments:
- `plan`:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D!(ccf_out::A1, ccf_covar_out::A2, λ::A3, flux::A4, var::A5,
                plan::PlanT = BasicCCFPlan(); projection_workspace::AbstractArray{T6,3} = zeros(length(λ),1,size(ccf_out,1)),
                ccf_var_scale::Real = 1.0, assume_sorted::Bool = false, calc_covar_terms::Bool = false ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,1}, T4<:Real, A4<:AbstractArray{T4,1}, T5<:Real, A5<:AbstractArray{T5,1}, T6<:Real,
                PlanT<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    @assert assume_sorted || issorted( plan.line_list.λ )
    v_grid = calc_ccf_v_grid(plan)
    num_vels = length(v_grid)
    @assert size(ccf_out,1) == num_vels
    @assert size(ccf_covar_out,1) == size(ccf_covar_out,2)
    @assert size(ccf_covar_out,1) == num_vels


    if length(plan.line_list) < 1  # If no lines in chunk
        ccf_out .= zero(eltype(A1))
        ccf_covar_out .= zero(eltype(A2))
        return (ccf=ccf_out, ccf_var=ccf_covar_out)
    end

    @assert size(projection_workspace,1) == length(λ)
    @assert size(projection_workspace,2) == 1   # Keep until implement CCF w/ multiple masks simultaneously
    @assert size(projection_workspace,3) == size(ccf_out,1)

    ccf_covar_out .= zero(eltype(A2))
    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:num_vels
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(view(projection_workspace, :, :, i), λ, plan, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = plan.allow_nans ? NaNMath.sum(flux .* view(projection_workspace,:,1,i) ) : sum(flux .* view(projection_workspace,:,1,i) )
        # computing the variance term of CCF covar's diagonal entries moved to separate loop for memory access locality
        # ccf_covar_out[i,i] = plan.allow_nans ? NaNMath.sum(var .* view(projection_workspace,:,1,i).^2) : sum(var .* view(projection_workspace,:,1,i).^2)
    end
    # Add variance terms for off-diagonal terms of CCF covar
    for i in 1:num_vels
        for j in i:num_vels
            if i == j
                ccf_covar_out[i,i] = plan.allow_nans ? NaNMath.sum(var .* view(projection_workspace,:,1,i).^2) : sum(var .* view(projection_workspace,:,1,i).^2)
            else  # i != j
                covar_term = plan.allow_nans ? NaNMath.sum(var .* view(projection_workspace,:,1,i) .* view(projection_workspace,:,1,j) ) : sum(var .* view(projection_workspace,:,1,i) .* view(projection_workspace,:,1,j) )
                ccf_covar_out[i,j] = covar_term
                ccf_covar_out[j,i] = covar_term
            end  # if i== j ... else ...
        end # for j
    end # for i
    #= =#
    if calc_covar_terms   # If we want to include covar terms due to LSF  # WARNING: SLOW
        pdf_norm(x,σ²) = exp(-0.5*x^2/σ²)/sqrt(2π*σ²)
        lsf_σ = 2400.0
        ccf_off_diag_σ_factor = 2.0

        for k in 1:length(λ)
            for l in (k+1):length(λ)
                #if k == l    continue   end     #  l starts at k+1 anyway; variance terms computed above
                Δv = log(λ[k]/λ[l])*speed_of_light_mps
                if abs(Δv) > (ccf_off_diag_σ_factor * lsf_σ)   continue    end
                flux_covar_norm = pdf_norm(Δv,2*lsf_σ^2)
                ave_var_λ = sqrt(var[k] * var[l])    # TODO: Improve this approximation based on flux var at v_mean?
                prefactor = flux_covar_norm * ave_var_λ
                for i in 1:num_vels
                    for j in i:num_vels
                        covar_term = prefactor * projection_workspace[k,1,i] * projection_workspace[l,1,j]
                        if i == j         # diagnonal terms of covar
                            ccf_covar_out[i,i] += 2*covar_term
                        else # i != j     # off diagonal terms of covar
                            ccf_covar_out[i,j] += 2*covar_term
                            ccf_covar_out[j,i] += 2*covar_term
                        end   # if i==j ... else ...
                    end # for j
                end # for i
            end # for l
        end  # for k
    else  #  !calc_covar_terms
        ccf_covar_out .*= ccf_var_scale
    end # if false
    return (ccf=ccf_out, ccf_covar=ccf_covar_out )
end


"""
    `ccf_1D( λs, fluxes; ccf_plan )`
    Compute the cross correlation function of a spectrum with a mask.
# Inputs:
- `λs`: 1-d array of wavelengths
- `fluxes`:  1-d array of fluxes
# Optional Arguments:
- `ccf_plan`:  parameters for computing ccf (BasicCCFPlan())
# Returns:
- 1-d array of size length(calc_ccf_v_grid(plan))
"""
function ccf_1D(λ::A2, flux::A3,
                plan::PlanT = BasicCCFPlan() ) where {
                T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                PlanT<:AbstractCCFPlan }
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    len_v_grid = calc_length_ccf_v_grid(plan)
    ccf_out = zeros(len_v_grid)
    if length(plan.line_list) < 1  # If no lines in chunk
        return ccf_out
    end
    ccf_1D!(ccf_out, λ, flux, plan )
    return ccf_out
end

"""
    `ccf_1D( λs, fluxes; ccf_plan )`
    Compute the cross correlation functions of spectra with a mask.
# Inputs:
- `λs`: 1-d array of wavelengths
- `fluxes`:  2-d array of fluxes, individual spectra along first dim
# Optional Arguments:
- `ccf_plan`:  parameters for computing ccf (BasicCCFPlan())
# Returns:
- 2-d array of size (length(calc_ccf_v_grid(plan)), size(flux, 2))
"""
function ccf_1D(λ::A2, flux::A3,
                plan::PlanT = BasicCCFPlan() ) where {
                T2<:Real, A2<:AbstractArray{T2,1},
                T3<:Real, A3<:AbstractArray{T3,2},
                PlanT<:AbstractCCFPlan }
    @assert ndims(λ) == 1
    @assert ndims(flux) == 2
    @assert length(λ) == size(flux, 1)
    len_v_grid = calc_length_ccf_v_grid(plan)
    ccf_out = zeros(len_v_grid, size(flux, 2))
    if length(plan.line_list) < 1  # If no lines in chunk
        return ccf_out
    end

    for i in 1:size(flux, 2)
        ccf_1D!(view(ccf_out, :, i), λ, view(flux, :, i), plan)
    end
    return ccf_out
end

"""
    `ccf_1D( λs, fluxes, vars; ccf_plan )`

Compute the cross correlation function of a spectrum with a mask and its variance^1.
^1 = This version computes the diagonal terms for the ccf variance and neglects covariances due to the line spread function.
Can roughly compensate by scaling up the ccf_var_scale from the default of 1.
WARNING: Still needs more testing.
# Inputs:
- `λs`: 1-d array of wavelengths
- `fluxes`:  1-d array of fluxes
# Optional Arguments:
- `ccf_plan`:  parameters for computing ccf (BasicCCFPlan())
# Returns:
- 1-d array of size length(calc_ccf_v_grid(plan))
"""
function ccf_1D(λ::A2, flux::A3, var::A4,
                plan::PlanT = BasicCCFPlan(); ccf_var_scale::Real = 1.0,
                calc_ccf_covar::Bool = false, calc_covar_terms::Bool = false ) where {
                T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1}, T4<:Real, A4<:AbstractArray{T4,1},
                PlanT<:AbstractCCFPlan }
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    len_v_grid = calc_length_ccf_v_grid(plan)
    ccf_out = zeros(len_v_grid)
    ccf_covar_out = calc_ccf_covar ? zeros(len_v_grid,len_v_grid) : zeros(len_v_grid)
    if length(plan.line_list) < 1  # If no lines in chunk
        return (ccf=ccf_out, ccf_var=ccf_covar_out)
    end
    if calc_ccf_covar
        ccf_1D!(ccf_out, ccf_covar_out, λ, flux, var, plan, ccf_var_scale=ccf_var_scale, calc_covar_terms=calc_covar_terms )
        return (ccf=ccf_out, ccf_covar=ccf_covar_out)
    else
        ccf_1D!(ccf_out, ccf_covar_out, λ, flux, var, plan, ccf_var_scale=ccf_var_scale )
        return (ccf=ccf_out, ccf_var=ccf_covar_out)
    end
end


""" `project_mask!( output, λs, ccf_plan; shift_factor )`

Compute the projection of the mask onto the 1D array of wavelengths (λs) at a given shift factor (default: 1).
The mask is computed from the ccf_plan, including a line_list and mask_shape (default: tophat).
Assumes plan.line_list is already sorted.

"""
function project_mask!(projection::A2, λs::A1, plan::PlanT ; shift_factor::Real=one(T1)
                            ) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, PlanT<:AbstractCCFPlan }
    num_lines_in_mask = length(plan.line_list)
    @assert num_lines_in_mask >= 1
    @assert size(projection,1) == length(λs)
    @assert size(projection,2) == 1         # Keep until implement CCF w/ multiple masks simultaneously
    projection .= zero(eltype(projection))  # Important to zero out projection before accumulating there!

    # set indexing variables
    p = 1
    p_left_edge_of_current_line = 1
    λsle_cur = λs[p]-0.5*(λs[p+1]-λs[p])   # Current pixel's left edge
    λsre_cur = 0.5*(λs[p]+λs[p+1])   # Current pixel's left edge
    #@assert issorted( plan.line_list.λ )  # Moved outside of inner loop
    m = searchsortedfirst(plan.line_list.λ, λsle_cur/shift_factor)
    while m <= num_lines_in_mask && λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor < λsle_cur   # only includemask entry if it fit entirely in chunk
        m += 1
    end
    if m > num_lines_in_mask   # Return early if no lines fall within chunk's lambda's
        return projection
    else
        #println("# Starting with mask entry at λ=", plan.line_list.λ[m], " for chunk with λ= ",first(λs)," - ",last(λs))
    end
    on_mask = false
    mask_mid = plan.line_list.λ[m] * shift_factor
    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's upper limit
    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's lower limit

    c_mps = speed_of_light_mps
    mask_weight = plan.line_list.weight[m]
    # loop through lines in mask, weight mask by fraction of PSF falling in each wavelength bin and divide by Δlnλ for pixel
    while m <= num_lines_in_mask
        if !on_mask
            if λsre_cur > mask_lo
                if λsre_cur > mask_hi   # Pixel overshot this mask entry, so weight based on mask filling only a portion of the pixel,
                    #For Tophat: projection[p,1] += mask_weight * (mask_hi - mask_lo) / (λsre_cur - λsle_cur)
                    #frac_of_psf_in_v_pixel = 1.0
                    one_over_delta_z_pixel = 0.5*(λsre_cur + λsle_cur) / (λsre_cur - λsle_cur)
                    projection[p,1] += one_over_delta_z_pixel * mask_weight
                    m += 1
                    if m<=num_lines_in_mask      # Move to next line
                        mask_mid = plan.line_list.λ[m] * shift_factor
                        mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_weight = plan.line_list.weight[m]
                    else         # We're out of lines, so exit
                        break
                    end
                else                       # Right edge of current pixel has entered the mask for this line, but left edge hasn't
                    #For Tophat: projection[p,1] += mask_weight * (λsre_cur - mask_lo) / (λsre_cur - λsle_cur)
                    frac_of_psf_in_v_pixel = integrate(plan.mask_shape, (mask_lo-mask_mid)/mask_mid*c_mps, (λsre_cur-mask_mid)/mask_mid*c_mps)
                    one_over_delta_z_pixel = 0.5*(λsre_cur + λsle_cur) / (λsre_cur - λsle_cur)
                    projection[p,1] += frac_of_psf_in_v_pixel * one_over_delta_z_pixel * mask_weight
                    on_mask = true         # Indicate next pixel will hav something to contribute based on the current line
                    p_left_edge_of_current_line = p                 # Mark the starting pixel of the line in case the lines overlap
                    p+=1                   # Move to next pixel
                    λsle_cur = λsre_cur   # Current pixel's left edge
                    λsre_cur = p<size(projection,1) ? 0.5*(λs[p]+λs[p+1]) : λs[p]+0.5*(λs[p]-λs[p-1])    # Current pixel's right edge

                end
            else   # ALl of pixel is still to left of beginning of mask for this line.
                p+=1    # TODO: Opt try to accelerate with something like searchsortedfirst?  Probably not worth it.
                λsle_cur = λsre_cur   # Current pixel's left edge
                λsre_cur = p<size(projection,1) ? 0.5*(λs[p]+λs[p+1]) : λs[p]+0.5*(λs[p]-λs[p-1])    # Current pixel's right edge
            end
        else  # on_mask == true
            if λsre_cur > mask_hi               # Right edge of this pixel moved past the edge of the current line's mask
                #For Tophat: projection[p,1] += mask_weight * (mask_hi - λsle_cur) / (λsre_cur - λsle_cur)
                frac_of_psf_in_v_pixel = integrate(plan.mask_shape, (λsle_cur-mask_mid)*c_mps/mask_mid, (mask_hi-mask_mid)*c_mps/mask_mid)
                one_over_delta_z_pixel = 0.5*(λsre_cur + λsle_cur) / (λsre_cur - λsle_cur)
                projection[p,1] += frac_of_psf_in_v_pixel * one_over_delta_z_pixel * mask_weight
                on_mask = false                 # Indicate that we're done with this line
                m += 1                          # Move to next line
                if m<=num_lines_in_mask
                    mask_mid = plan.line_list.λ[m] * shift_factor
                    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_weight = plan.line_list.weight[m]
                    if mask_lo < λsle_cur       # If the lines overlap, we may have moved past the left edge of the new line. In that case go back to the left edge of the previous line.
                        p = p_left_edge_of_current_line
                        λsle_cur = p>1 ? 0.5*(λs[p]+λs[p-1]) : λs[p]-0.5*(λs[p+1]-λs[p])   # Current pixel's left edge
                        λsre_cur = 0.5*(λs[p]+λs[p+1])   # Current pixel's right edge
                    end

                else                            # We're done with all lines, can return early
                    break
                end
            else                                # This pixel is entirely within the mask
                # assumes plan.mask_shape is normalized to integrate to unity and flux is constant within pixel
                #For Tophat: projection[p,1] += mask_weight    # Add the full mask weight
                frac_of_psf_in_v_pixel = integrate(plan.mask_shape, (λsle_cur-mask_mid)/mask_mid*c_mps, (λsre_cur-mask_mid)/mask_mid*c_mps)
                one_over_delta_z_pixel = 0.5*(λsre_cur + λsle_cur) / (λsre_cur - λsle_cur)
                projection[p,1] += frac_of_psf_in_v_pixel * one_over_delta_z_pixel * mask_weight
                p += 1                          # move to next pixel
                λsle_cur = λsre_cur   # Current pixel's left edge
                λsre_cur = p<size(projection,1) ? 0.5*(λs[p]+λs[p+1]) : λs[p]+0.5*(λs[p]-λs[p-1])    # Current pixel's right edge
            end
        end
        if p>=size(projection,1) break end       # We're done with all pixels in this chunk.
    end
    projection ./= c_mps
    return projection
end
