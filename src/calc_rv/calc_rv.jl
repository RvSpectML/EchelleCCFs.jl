"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored and optimized by Eric Ford
"""

"""  Module for estimating the radial velocity based on the CCF  """
module RVFromCCF

import Statistics: mean
import PDMats: PDiagMat
#import Polynomials
import ..RvSpectMLBase: searchsortednearest, speed_of_light_mps


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

import ..EchelleCCFs: calc_ccf_template, calc_normalized_ccfs
include("ccf_sample_covar.jl")
export calc_ccf_sample_covar

include("ccf_covar_model.jl")
export calc_ccf_sample_covar_reduced_rank, calc_ccf_sample_covar_and_near_diag_model, est_covar_for_obs

include("fit_template_plus_deriv_to_ccf.jl")
#import .FitTemplateToCCF
#import .FitTemplateToCCF: MeasureRvFromCCFTemplate
export MeasureRvFromCCFTemplate

include("fit_template_plus_deriv_to_ccf_non_diag_covar.jl")
#import .FitTemplateToCCF
#import .FitTemplateToCCF: MeasureRvFromCCFTemplateNonDiagCovar
export MeasureRvFromCCFTemplateNonDiagCovar

end # module
