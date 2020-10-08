#using NamedTupleTools

function add_tuple_sum(x::NT, y::NT) where { NT<:NamedTuple }
    @assert all(fieldtypes(x) .== fieldtypes(y))
    names_of_values = propertynames(x)
    @assert propertynames(y)  == names_of_values
    values_of_x = fieldvalues(x)
    values_of_y = fieldvalues(y)
    #nt = namedtuple(names_of_values, values_of_x .+ values_of_y )
    values = map(t->t[1].+t[2], zip(values_of_x,values_of_y))
    nt = namedtuple(names_of_values, values )
    return nt
end
