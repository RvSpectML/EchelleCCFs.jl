"""
    Code to compute CCF
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored, optimized and generalized by Eric Ford
"""

""" Module for computing CCFs """
module EchelleCCFs

import RvSpectMLBase
import RvSpectMLBase: AbstractChunkOfSpectrum, AbstractChunkList, AbstractChunkListTimeseries
import RvSpectMLBase: AbstractInstrument
import RvSpectMLBase: num_chunks

#import Polynomials

include("physics.jl")

# Default values shared across instruments
const default_line_width_mps = RvSpectMLBase.default_line_width_mps  # m/s
default_chunk_size_factor = 3        # For default_calc_chunk_width TODO: Figure out what value to use.  Ask Alex
default_min_chunk_Δv = 20000           # m/s  for ChunkWidthFixedΔlnλ

default_Δv_to_avoid_tellurics = 30000  # m/s

include("masks/io.jl")
export AbstractLineList, BasicLineList

include("masks/masks.jl")

include("mask_shapes/mask_shapes.jl")
# exports its own types

# default parameters for CCF plans
const default_v_center = 0.0    # m/s
const default_v_step = 250.0    # m/s
const default_v_max = 15.0e3    # m/s
const default_v_width = 410.0   # m/s
const default_v_range_no_mask_change = default_v_max # m/s

include("ccf/ccf.jl")
# exports its own types/functions

include("calc_rv/calc_rv.jl")
using .RVFromCCF
export measure_rv_from_ccf, measure_rvs_from_ccf
export MeasureRvFromMinCCF, MeasureRvFromCCFCentroid, MeasureRvFromCCFQuadratic, MeasureRvFromCCFGaussian

end
