"""
    Code to loading mask shapes for CCF
Author: Eric Ford
Created: September 2020
"""

# constants

default_gaussian_ccf_σ = 5000.0
default_gaussian_ccf_truncation_scale_factor = 2*(sqrt(2*log(2)))
default_supergaussian_ccf_fwhm = 448.0*2*sqrt(2*log(2)) #*4.5
default_supergaussian_ccf_exponent = 1.3
default_supergaussian_ccf_truncation_scale_factor = 2
default_gaussian_mixture_ccf_truncation_Δv = 20000.0

"A struct implementing a specific mask shapes should be a subtype of AbstractCCFMaskShape."
abstract type AbstractCCFMaskShape  end
export AbstractCCFMaskShape

include("tophat.jl")
export TopHatCCFMask
include("gaussian.jl")
export GaussianCCFMask
include("halfcos.jl")
export CosCCFMask
include("supergaussian.jl")
export SuperGaussianCCFMask

# TODO: Fix before adding
#include("gaussian_mixture.jl")
#export GaussianMixtureCCFMask
