"""
Code for creating, loading and maniuplating line lists and masks.

Author: Eric Ford
Created: August 2020
"""

"""
    Read line list in ESPRESSO csv format.
ESPRESSO format: lambda and weight.
convert_air_to_vacuum determines whether to convert to vacuum wavelengths.
Warning: ESPRESSO masks don't provide line depth and sometimes include one entry for a blend of lines.
"""
function read_linelist_espresso(fn::String; convert_air_to_vacuum::Bool = true)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","weight"],delim=' ',ignorerepeated=true)
    @assert hasproperty(df, :lambda)
    @assert hasproperty(df, :weight)
    if convert_air_to_vacuum
        df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda])
    end
    return df
end

""" Read line list in VALD csv format.
   VALD format: lambda_lo, lambdaa_hi and depth.
   convert_air_to_vacuum determines whether to convert to vacuum wavelengths.
"""
function read_linelist_vald(fn::String; convert_air_to_vacuum::Bool = true)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","depth"])
    @assert hasproperty(df, :lambda)
    @assert hasproperty(df, :depth)
    if convert_air_to_vacuum
        df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda])
    end
    df[!,:weight] = df[!,:depth] # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    return df
end

""" Read line list in csv format.
   format: lambda, weight, lambda_lo, lambdaa_hi.
   Assumes air to vacuumb wavelength conversion has already been applied.
"""
function read_linelist_rvspectml(fn::String)
    local df = CSV.read(fn,DataFrame,threaded=false)
    @assert hasproperty(df, :lambda) || hasproperty(df, :wavelength)
    @assert hasproperty(df, :weight) || hasproperty(df, :depth) 
    if !hasproperty(df, :lambda) && hasproperty(df, :wavelength)
      rename!(df,:wavelength => :lambda)
    end
    if hasproperty(df, :depth) && !hasproperty(df, :weight)
        df[!,:weight] = df[!,:depth] # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    end
    return df
end
