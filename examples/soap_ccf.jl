make_plots = true

using EchelleCCFs
using Statistics
using DataFrames, Query

include("read_soap_spectra.jl")

λ_min = max(3950, minimum(λ) )  # Avoid Ca H & K lines
λ_max = min(6800, maximum(λ) )  # Avoid worst of telluric regions
idx_min = findfirst(x->x>=λ_min, λ)
idx_max = findlast(x->x<=λ_max, λ)
λ = λ[idx_min:idx_max]
flux = flux[idx_min:idx_max]

flux ./= mean(flux)

if make_plots
    using Plots
    plt = plot()
    #plt_idx = 1:length(λ)
    plt_idx = 1:2000
    plot!(plt, λ[plt_idx], flux[plt_idx], label=:none)
end


line_list_filename = joinpath(pkgdir(EchelleCCFs),"data","masks","G2.espresso.mas")
espresso_df = EchelleCCFs.read_linelist_espresso(line_list_filename)# # |>
                #@filter(λ_min < _.lambda < λ_max) |> DataFrame
line_list = EchelleCCFs.BasicLineList(espresso_df.lambda, espresso_df.weight)

resolution = 137000  # approximate R for EXPRESS
mask_width = EchelleCCFs.speed_of_light_mps / resolution  # m/s  Using a Gaussian CCF mask approximates effects of instrument resolution (neglects finite resolution of FTS)
mask_shape = GaussianCCFMask(mask_width)
#mask_shape = TopHatCCFMask(mask_width)   # In you wanted tophat mask for some reason
approx_v_offset_espresso_to_soap = -4.5e3  # just by eye so things aren't wildly shifted
max_bc = 30e3  #  ~ 2pi AU/year  in m/s
Δv_step = 400  # m/s arbitrary probably smaller/slower than you need
Δv_max = 30e3  # m/s arbitrary
ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max, midpoint=approx_v_offset_espresso_to_soap,
    range_no_mask_change=max_bc)
v_grid = calc_ccf_v_grid(ccf_plan) .- approx_v_offset_espresso_to_soap
ccf = ccf_1D( λ, flux, ccf_plan)

if make_plots
    using Plots
    plot(v_grid,ccf)
    scatter(v_grid,ccf,markersize=1.2,label="")
    xlabel!("Δv (m/s)")
    ylabel!("CCF ")
end
