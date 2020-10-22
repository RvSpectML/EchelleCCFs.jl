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
import RvSpectMLBase: AbstractSpectra, AbstractSpectra1D, AbstractSpectra2D
import RvSpectMLBase: AbstractChunkOfSpectrum, AbstractChunkList, AbstractChunkListTimeseries
import RvSpectMLBase: AbstractInstrument
import RvSpectMLBase: num_chunks

using DataFrames, Query, CSV

#import Polynomials

include("physics.jl")

# Default values shared across instruments
const default_line_width_mps = RvSpectMLBase.default_line_width_mps  # m/s
default_chunk_size_factor = 3        # For default_calc_chunk_width TODO: Figure out what value to use.  Ask Alex
default_min_chunk_Δv = 20000           # m/s  for ChunkWidthFixedΔlnλ

default_Δv_to_avoid_tellurics = 30000  # 2pi AU /year in m/s with a a little extra padding for Earth's rotation

include("masks/io.jl")
export AbstractLineList, BasicLineList

include("masks/masks.jl")

include("mask_shapes/mask_shapes.jl")
# exports its own types

include("line_list/line_list.jl")
export AbstractLineList, BasicLineList
export assign_lines_to_orders, calc_snr_weights_for_lines!

# default parameters for CCF plans
const default_v_center = 0.0    # m/s
const default_v_step = 250.0    # m/s
const default_v_max = 15.0e3    # m/s
const default_v_width = 410.0   # m/s
const default_v_range_no_mask_change = default_v_max # m/s

include("ccf/ccf.jl")
# exports its own types/functions

include("convenience/convenience.jl")
export calc_ccf_chunk, calc_ccf_chunk!
export calc_ccf_and_var_chunk, calc_ccf_and_var_chunk!
export calc_ccf_chunklist #, calc_ccf_chunklist!
export calc_ccf_and_var_chunklist #, calc_ccf_and_var_chunklist!
export calc_order_ccfs_chunklist, calc_order_ccfs_chunklist!
export calc_order_ccf_and_vars_chunklist, calc_order_ccf_and_vars_chunklist!
export calc_ccf_chunklist_timeseries
export calc_ccf_and_var_chunklist_timeseries

include("calc_rv/calc_rv.jl")
using .RVFromCCF
export RVFromCCF
export measure_rv_from_ccf, measure_rvs_from_ccf
export MeasureRvFromMinCCF, MeasureRvFromCCFCentroid, MeasureRvFromCCFQuadratic, MeasureRvFromCCFGaussian, MeasureRvFromCCFTemplate

end
