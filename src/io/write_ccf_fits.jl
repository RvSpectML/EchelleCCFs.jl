#using FITSIO
#using Pkg
#using Dates

function make_ccf_fits_header(metadata::AbstractDict; line_list_filename::Union{Nothing,String} = "unknown",
            v::Union{Nothing,Float64} = nothing, e_v::Union{Nothing,Float64} = nothing,
            blue_ord::Union{Nothing,Int64} = nothing, red_ord::Union{Nothing,Int64} = nothing,
            mask_shape::Union{Nothing,String} = nothing, wave_cal::Union{Nothing,String} = nothing,
            blazenorm::Union{Nothing,Bool} = nothing, contnorm::Union{Nothing,Bool} = nothing,
            div_tell::Union{Nothing,Bool} = nothing, rawnorm::Union{Nothing,Bool} = nothing
            )
    header_dict = Dict{String,Tuple{Any,String}}("VERSION"=>(string(Pkg.project().version), string(Pkg.project().name) * " code version"))
    # Add version info for packages used
    deps = Pkg.dependencies()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        dep.version === nothing && continue
        header_dict["VERSION_" * dep.name] = (string(dep.version), "version of package dependancy")
    end
    # header fields similar to EXPRES CCFs
    header_dict["DATE-CCF"] = (string(Dates.now()), "Date and time of CCF calculation")
    if hasproperty(metadata,:Filename)
        header_dict["FILENAME"] = (last(splitpath(metadata[:Filename])), "Original filename")
    end
    if hasproperty(metadata,:bjd)
        header_dict["MJD"] = (metadata[:bjd], "MJD of geometric midpoint")
    end
    if !isnothing(line_list_filename)
        header_dict["MASK"] = (line_list_filename, "CCF mask linelist used")
    end
    # Skipped: EXPCOUNT, SNR, V, E_V, CHI2, BLUE_ORD, RED_ORD
    if !isnothing(mask_shape)
        header_dict["WINDOW"] = (mask_shape, "Shape of CCF mask used in CCF")
    end
    # Skipped BARYCOR
    if hasproperty(metadata,:airmass)
        header_dict["AIRMASS"] = (metadata[:airmass], "Airmass of exposure")
    end
    # Skipped EXPTIME
    if hasproperty(metadata,:expres_epoch)
        header_dict["RVEPOCH"] = (metadata[:expres_epoch], "RV calibration epoch")
    end
    if hasproperty(metadata,:wave_cal)
        header_dict["WAVE_CAL"] = (wave_cal, "Wavelength calibration used when computing CCF")
    end
    # Skipped CCOR
    if !isnothing(blazenorm)
        header_dict["BLAZE"] = (blazenorm, "Pixels weighted by blaze in CCF")
    end
    if !isnothing(contnorm)
        header_dict["CONTNORM"] = (contnorm, "Continuum normalized spectrum in CCF")
    end
    # Skipped SUB_CONT
    if !isnothing(div_tell)
        header_dict["DIV_TELL"] = (div_tell, "Tellurics divided before CCF")
    end
    # Skipped PLATESCL
    # Added some that I thought might be useful to have
    if hasproperty(metadata,:normalization)
        header_dict["NORMALIZATION"] = (string(metadata[:normalization]), "How spectrum was normalized before computing CCF")
    end
    if hasproperty(metadata,:snr_prelim)
        header_dict["SNR_PRELIM"] = (string(metadata[:snr_prelim]), "Something related to SNR, but different from what EXPRES provided")
    end
    if hasproperty(metadata,:pwv)
        header_dict["PWV"] = (metadata[:pwv], "precip. water vabor")
    end
    if hasproperty(metadata,:sundist)
        header_dict["SUN_DIST"] = (metadata[:sundist], "Distance to sun")
    end
    if hasproperty(metadata,:moondist)
        header_dict["MOON_DIST"] = (metadata[:moondist], "Distance to moon")
    end
    fits_header = FITSHeader(collect(keys(header_dict)),first.(values(header_dict)),last.(values(header_dict)))
end


function write_ccf_fits(filename::String, v_grid::AbstractArray{T1,1}, ccf_total::AbstractArray{T2,1}, ccf_total_var::AbstractArray{T3,1};
            orders::AbstractArray{T4,1}=zeros(Int64,0), ccf_orders::AbstractArray{T5,2}=zeros(T2,0,0), ccf_orders_var::AbstractArray{T6,2}=zeros(T3,0,0),
            fits_hdr::FITSHeader,
            #total_vel::Union{Nothing,Real} = nothing, sigma_total_vel::Union{Nothing,Real} = nothing,
            order_vels::Union{Nothing,AbstractArray{T7,1}} = nothing, sigma_order_vels::Union{Nothing,AbstractArray{T8,1}} = nothing ) where {
                T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real, T6<:Real, T7<:Real, T8<:Real   }
    f = FITS(filename, "w")
    #if isnothing(total_vel) || isnothing(sigma_total_vel)
    write(f,UInt8[],header=fits_hdr,name="metadata")
    #else
    #    write(f,UInt8[],header=fits_hdr,name="metadata",v=total_vel,e_V=sigma_total_vel)
    #end
    t2 = OrderedDict("V_grid"=>collect(v_grid),"ccf"=>ccf_total,"e_ccf"=>sqrt.(ccf_total_var) )
    write(f,collect(keys(t2)),collect(values(t2)),name="ccf_total")
    if length(orders) > 0 && length(ccf_orders) > 0 && length(ccf_orders_var) > 0
        if isnothing(order_vels) || isnothing(sigma_order_vels)
            t3 = OrderedDict("orders"=>convert.(Int64,orders),"ccfs"=>ccf_orders,"errs"=>sqrt.(ccf_orders_var) )
        else
            t3 = OrderedDict("orders"=>convert.(Int64,orders),"ccfs"=>ccf_orders,"errs"=>sqrt.(ccf_orders_var),
                                "v" => order_vels,"e_v" => sigma_order_vels)
        end
        write(f,collect(keys(t3)),collect(values(t3)),name="ccf_order")
    end
    close(f)
end

function write_each_ccf_fits(metadata::AbstractArray{Dict{Symbol,Any},1}, v_grid::AbstractArray{T1,1}, ccf_total::AbstractArray{T2,2}, ccf_total_var::AbstractArray{T3,2};
            orders::AbstractArray{T4,1}=zeros(Int64,0), ccf_orders::AbstractArray{T5,3}=zeros(T2,0,0,size(ccf_total,2)), ccf_orders_var::AbstractArray{T6,3}=zeros(T3,0,0,size(ccf_total,2)),
            fits_hdr::FITSHeader, #  )
            total_vel::Union{Nothing,AbstractArray{T7,1}} = nothing, sigma_total_vel::Union{Nothing,AbstractArray{T7,1}} = nothing,
            order_vels::Union{Nothing,AbstractArray{T8,2}} = nothing, sigma_order_vels::Union{Nothing,AbstractArray{T8,2}} = nothing ) where {
                T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real, T6<:Real, T7<:Real, T8<:Real }
    @assert length(metadata) == size(ccf_total,2) == size(ccf_total_var,2)
    num_files = length(metadata)
    for i in 1:num_files
        fits_fn = last(splitpath(metadata[i][:Filename]))
        ccf_fn = replace(fits_fn, ".fits" => "_ccf.fits")
        if isnothing(total_vel) || isnothing(sigma_total_vel)
            fits_hdr = make_ccf_fits_header(metadata[i])
        else
            fits_hdr = make_ccf_fits_header(metadata[i],v=total_vel[i],e_v=sigma_total_vel[i])
        end
        if isnothing(order_vels) || isnothing(sigma_order_vels)
            write_ccf_fits(ccf_fn, v_grid, ccf_total[:,i], ccf_total_var[:,i],
                orders=orders, ccf_orders=ccf_orders[:,:,i], ccf_orders_var=ccf_orders_var[:,:,i],
                fits_hdr=fits_hdr)
        else
            @assert size(order_vels) == size(sigma_order_vels)
            @assert size(order_vels,1) == num_files
            @assert size(order_vels,2) == length(orders)
            write_ccf_fits(ccf_fn, v_grid, ccf_total[:,i], ccf_total_var[:,i],
                orders=orders, ccf_orders=ccf_orders[:,:,i], ccf_orders_var=ccf_orders_var[:,:,i],
                fits_hdr=fits_hdr,
                order_vels=order_vels[i,:], sigma_order_vels=sigma_order_vels[i,:] )
        end
        # break  # WARNING: Aborts after one file for testing purposes
    end
end
