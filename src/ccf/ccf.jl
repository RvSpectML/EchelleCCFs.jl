# default parameters
const default_v_center = 0.0    # m/s
const default_v_step = 250.0    # m/s
const default_v_max = 15.0e3    # m/s
const default_v_width = 410.0   # m/s
const default_v_range_no_mask_change = default_v_max # m/s


include("plan.jl")
export AbstractCCFPlan, BasicCCFPlan
export calc_ccf_v_grid, calc_length_ccf_v_grid

include("calc_ccf.jl")
export ccf_1D, ccf_1D!

include("ccf_template.jl")
export calc_normalized_ccfs, calc_ccf_template

#=
include("convenience/convenience.jl")
export calc_ccf_chunk, calc_ccf_chunk!
export calc_ccf_and_var_chunk, calc_ccf_and_var_chunk!
export calc_ccf_chunklist #, calc_ccf_chunklist!
export calc_ccf_and_var_chunklist #, calc_ccf_and_var_chunklist!
export calc_order_ccfs_chunklist, calc_order_ccfs_chunklist!
export calc_order_ccf_and_vars_chunklist, calc_order_ccf_and_vars_chunklist!
export calc_ccf_chunklist_timeseries
export calc_ccf_and_var_chunklist_timeseries
=#

# Will soon comment out to reduce precompilation time
#include("calc_ccf_old.jl")
#include("convenience_old.jl")
#export ccf_1D_old, ccf_1D_old!
