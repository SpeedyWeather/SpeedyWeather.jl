export tree

"""
$(TYPEDSIGNATURES)
Create a tree of fields inside a Simulation instance and fields within these fields
as long as they are defined within the modules argument (default SpeedyWeather).
Other keyword arguments are `max_level::Integer=10`, `with_types::Bool=false`."""
function tree(S::Simulation{M}; modules=SpeedyWeather, kwargs...) where M
    println("Simulation{$(model_type(M))}")
    _tree(S, modules; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Create a tree of fields inside a model and fields within these fields
as long as they are defined within the modules argument (default SpeedyWeather).
Other keyword arguments are `max_level::Integer=10`, `with_types::Bool=false`."""
function tree(M::ModelSetup; modules=SpeedyWeather, kwargs...)
    println("$(model_type(M)){...}")
    _tree(M, modules; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Create a tree of fields inside S and fields within these fields
as long as they are defined within the modules argument (default SpeedyWeather).
Other keyword arguments are `max_level::Integer=10`, `with_types::Bool=false`."""
function tree(S; modules=SpeedyWeather, kwargs...)
    println("$(typeof(S))")
    _tree(S, modules; kwargs...)
end

function _tree(
    S,
    modules::Module...;         # fields within the modules are further inspected
    max_level::Integer = 10,    # depth of the tree
    with_types::Bool = false,   # print also the types of the fields?
)
    level = 0                   # starting level of depth tree
    prevs = falses(max_level+2) # determine whether there's still a branch levels up
                                # needed for │ printing of upper levels that continue
    property_names = propertynames(S)
    n_properties = length(property_names)

    for (i,property_name) in enumerate(property_names)
        last = i == n_properties    # last elements in branches are printed with └ not ├
        prevs[level+1] = ~last
        print_branch(property_name, S, level, max_level, last, prevs, with_types, modules...)
    end
end

function print_branch(
    property_name::Symbol,  # name of current field
    S,                      # its parent struct
    level::Integer,         # depth of tree we're on
    max_level::Integer,     # maximum depth of tree
    last::Bool,             # is that the last field in parent?
    prevs::BitVector,       # branching of previous fields completed?
    with_types::Bool,       # print also types
    modules::Module...,     # in which modules to search
 )
    level == max_level && return nothing

    property = getproperty(S,property_name)
    child = parentmodule(typeof(property)) in modules

    continue_branching = false
    is_array = false
    if child
        continue_branching = true
    elseif property isa AbstractArray
        if length(property) > 0
            property = property[1]      # unpack first element in array
            continue_branching = parentmodule(typeof(property)) in modules
            is_array = continue_branching
        end
    end

    property_names = propertynames(property)
    n_properties = length(property_names)
    junction1 = (n_properties > 0 && child) || is_array ? "┐" : " "

    vertical_lines = prod([b ? "│" : " " for b in prevs[1:level]])
    junction2 = last ? "└" : "├"
    print(vertical_lines, junction2, junction1, property_name)
    s = ~with_types ? "" : "::$(typeof(property))"
    println(s)

    if continue_branching
        for (i,branch) in enumerate(property_names)
            last = i == n_properties
            prevs[level+2] = ~last
            print_branch(branch, property, level+1, max_level, last, prevs, with_types, modules...)
        end
    end
end