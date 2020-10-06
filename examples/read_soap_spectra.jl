import Pkg
Pkg.activate(joinpath(pkgdir(EchelleCCFs),"examples"))

using HDF5, FileIO

include("../deps/download_soap_example.jl")

function read_soap_output_gilbertson_ford_dumuseque_2020(filename::String, obs_idx::Integer)
    local time, λ, flux
    c = h5open(filename, "r") do file
        flux_size = size(file["active"])
        @assert 1<= obs_idx <= flux_size[2]
        flux = file["active"][1:flux_size[1],obs_idx]
        λ =  read(file,"lambdas")
        idx_min = findfirst(x->x>0, flux)
        idx_max = findlast(x->x>0, flux)
        λ = view(λ,idx_min:idx_max)
        flux = view(flux,idx_min:idx_max)
        p_rot_eq = 25.05 # days
        time = file["phases"][obs_idx] * p_rot_eq
    end
    return (time=time, λ=λ, flux = flux)
end

filename = joinpath(pkgdir(EchelleCCFs),"data","spectra","res-1000-1years_full_id1.h5")
(time, λ, flux) = read_soap_output_gilbertson_ford_dumuseque_2020(filename, 1)
