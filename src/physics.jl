"""
Physics related to line lists and masks.

Author: Eric Ford
Created: August 2020
"""

# physical constants
const speed_of_light_mps = 2.99782458e8    # m/s


"""
    calc_doppler_factor(vel)

Compute the longitudinal relativistic doppler factor given a velocity
in meters per second.
"""
function calc_doppler_factor(vel::Real)
    one(vel) + vel/speed_of_light_mps
end

#=
 # Eric replaced this version with above.  Should double check that won't cause anyeone suprises.
function calc_doppler_factor(vel::Real)
    num = one(vel) + vel/c_ms
    den = one(vel) - vel/c_ms
    return sqrt(num/den)
end
=#

"""  Convert vacuum wavelength (in Å) to air wavelength
Ref: Donald Morton (2000, ApJ. Suppl., 130, 403) via
     https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
"""
function λ_vac_to_air(λ_vac::Real)
    @assert 3500 < λ_vac < 13000  # Making sure in Å for optical/NIR spectra.
    local s = 10000/λ_vac
    local n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    return λ_vac/n
end

""" Convert air wavelength (in Å) to vacuum wavelength
Ref: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
     VALD3 tools use the following solution derived by N. Piskunov
"""
function λ_air_to_vac(λ_air::Real)
    @assert 3500 < λ_air < 13000  # Making sure in Å for optical/NIR spectra.
    local s = 10000/λ_air
    local n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2)
    λ_air*n
end
