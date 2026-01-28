import Pkg
Pkg.activate(joinpath(pkgdir(EchelleCCFs),"examples"))

using EchelleCCFs
using HDF5, FileIO

function read_soap_demo(filename::String)
    local time, λ, flux
    c = h5open(filename, "r") do file
        λ =  read(file,"λ")
        flux =  read(file,"flux")
    end
    return (λ=λ, flux = flux)
end
filename = joinpath(pkgdir(EchelleCCFs),"data","spectra","soap_demo.h5")
if !isfile(filename)
    include(joinpath(pkgdir(EchelleCCFs), "deps", "build.jl")
end
(λ, flux) = read_soap_demo(filename)
