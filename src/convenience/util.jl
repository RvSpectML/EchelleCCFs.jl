#using NamedTupleTools

function add_tuple_sum(x::NT, y::NT) where { NT<:NamedTuple }
    @assert all(typeof(x).types .== typeof(y).types)
    names_of_values = keys(x)
    @assert keys(y)  == names_of_values
    values_of_x = values(x)
    values_of_y = values(y)
    #nt = namedtuple(names_of_values, values_of_x .+ values_of_y )
    values = map(t->t[1].+t[2], zip(values_of_x,values_of_y))
    nt = NamedTuple{names_of_values}( values )
    return nt
end
