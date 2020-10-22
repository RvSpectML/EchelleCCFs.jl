"""
Code for convenience functions for calculating CCFs on RvSpectMLBase's types
At some point could consider moving this to RvSpectML.

Author: Eric Ford
Created: August 2020
"""

using Statistics

using NamedTupleTools
include("util.jl")

include("ccf_chunk.jl")

include("ccf_chunklist.jl")

using ThreadedIterables
include("ccf_chunklist_timeseries.jl")

include("order_ccf_chunklist.jl")

#using ThreadTools
include("order_ccf_chunklist_timeseries.jl")
