"""
    Code for LineList types used by CCfs
Author: Eric Ford
Created: August 2020
"""

"""A struct implementing a line list should be a subtype of AbstractLineList. """
abstract type AbstractLineList end

""" A basic line list for passing to compute CCFs.
Contains (views into) arrays specifying the minimum and maximum wavelength range and weight for each line. """
struct BasicLineList{T<:Real, AA<:AbstractArray{T,1} } <: AbstractLineList
    λ::AA
    weight::AA
end

""" `BasicLineList( λ, weight )` """
function BasicLineList{T}(λ::AA, w::AA) where { T<:Real, AA<:AbstractArray{T,1} }
    @assert length(λ) == length(w)
    @assert length(λ) >= 1
    @assert all(0.0 .<= w .<= 1.0)
    BasicLineList{eltype(w),typeof(w)}(λ,w)
end

function BasicLineList( df::DataFrame )
    @assert hasproperty(df,:lambda)
    @assert hasproperty(df,:weight)
    BasicLineList{eltype(df[:,:lambda]),typeof(df[:,:weight])}(df[:,:lambda],df[:,:weight])
end

function EmptyBasicLineList()
    return BasicLineList{Float64,Array{Float64,1}}(zeros(0),zeros(0))
end


""" A basic line list for passing to compute CCFs.
Contains (views into) arrays specifying the minimum and maximum wavelength range and weight for each line. """
struct BasicLineList2D{T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Integer,  AA2<:AbstractArray{T2,1} } <: AbstractLineList
    λ::AA1
    weight::AA1
    order::AA2
end

""" `BasicLineList( λ, weight )` """
function BasicLineList2D{T1,T2}(λ::AA1, w::AA1, o::AA2) where {T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Integer,  AA2<:AbstractArray{T2,1} }
    @assert length(λ) == length(w) == length(o)
    @assert length(λ) >= 1
    @assert all(0.0 .<= w .<= 1.0)
    @assert all(1 .<= o .<= 200)           # arbitrary limit for maximum number of orders
    BasicLineList2D{eltype(w),typeof(w),eltype{o},typeof(o)}(λ,w,o)
end

function BasicLineList2D( df::DataFrame )
    @assert hasproperty(df,:lambda)
    @assert hasproperty(df,:weight)
    @assert hasproperty(df,:order)
    @assert all(0.0 .<= df[!,:weight] .<= 1.0)
    @assert all(1 .<= df[!,:order] .<= 200)           # arbitrary limit for maximum number of orders
    BasicLineList2D{eltype(df[:,:lambda]),typeof(df[:,:weight]),eltype(df[:,:order]),typeof(df[:,:order])}(df[!,:lambda],df[!,:weight],df[!,:order])
end

function EmptyBasicLineList2D()
    return BasicLineList2D{Float64,Array{Float64,1},Int32,Array{Int32,1}}(zeros(0),zeros(0), zeros(Int32,0) )
end


#=  Not fully implemented/tested yet
""" A line list for passing to compute CCFs with variable line widths.
Contains (views into) arrays specifying the minimum and maximum wavelength range and weight for each line. """
struct VarWidthLineList{T<:Real, AA<:AbstractArray{T,1} } <: AbstractLineList
    λ_lo::AA
    λ_hi::AA
    weight::AA
end

""" VarWidthLineList( λ_lo, λ_hi, weight ) """
function VarWidthLineList{T}(lo::AA, hi::AA, w::AA) where { T<:Real, AA<:AbstractArray{T,1} }
    @assert length(lo) == length(hi) == length(w)
    @assert length(lo) >= 1
    @assert 0.0 .<= λ_hi.-λ_lo .<= 2.0
    @assert 0.0 .<= w .<= 1.0
    VarWidthLineList{eltype(w),typeof(w)}(lo,hi,w)
end
=#

import Base.length
""" Return length of line list. """
length(ll::AbstractLineList) = length(ll.weight)
