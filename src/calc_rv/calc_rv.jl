"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored and optimized by Eric Ford
"""

"""
Module for estimating the radial velocity based on the CCF
"""
module RVFromCCF

import Statistics: mean
import PDMats: PDiagMat
#import Polynomials
import ..RvSpectMLBase: searchsortednearest

include("common.jl")
export AbstractMeasureRvFromCCF
export measure_rv_from_ccf, measure_rvs_from_ccf

include("find_min_ccf.jl")
export MeasureRvFromMinCCF

include("calc_ccf_centroid.jl")
export MeasureRvFromCCFCentroid

include("fit_quadratic_to_ccf.jl")
export MeasureRvFromCCFQuadratic

include("fit_gaussian_to_ccf.jl")
#import .FitGaussianToCCF
#import .FitGaussianToCCF: MeasureRvFromCCFGaussian
export MeasureRvFromCCFGaussian

end # module
