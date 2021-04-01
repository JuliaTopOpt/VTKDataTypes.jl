function similar_cells( cell_connectivity1::Vector{Int}, 
                        cell_type1::Int, 
                        cell_connectivity2::Union{Tuple{Vararg{Int}}, Vector{Int}},
                        cell_type2::Int
                      )
    for k in cell_connectivity2
        in(k, cell_connectivity1) || return false
    end
    cell_type1 == cell_type2 || return false
    return true
end

function coherent(a::T, b::T) where {T<:AbstractVTKMultiblockData}
    return same_ordered_geometry_shape(a,b) && same_data_shape(a,b)
end
function coherent(a::T, b::T) where {T<:AbstractVTKSimpleData}
    return same_geometry_shape(a, b) && same_data_shape(a, b)
end
function coherent(dataset::AbstractTimeSeriesVTKData)
    length(dataset) == 0 && return true 
    block1 = dataset[1]
    for block in dataset[2:end]
        coherent(block1, block) || return false
    end
    return true
end
function coherent(a::T, b::T) where {T<:AbstractTimeSeriesVTKData}
    length(a) == length(b) && 
    coherent(a) && coherent(b) && begin 
        _out = true
        for i in 1:length(a)
            (_out = coherent(a[i], b[i])) || break
        end
        _out
    end || return false
    return true
end
coherent(a::AbstractVTKData, b::AbstractVTKData) = false

extents(dataset::VTKStructuredData) = Base.tail(size(dataset.point_coords))
function extents(dataset::VTKRectilinearData)
    return ntuple(i->length(dataset.point_coords[i]), Val(length(dataset.point_coords)))
end
extents(dataset::VTKUniformRectilinearData) = dataset.extents
extents(dataset::AbstractVTKStructuredData, ind::Int) = extents(dataset)[ind]
extents(dataset::AbstractVTKStructuredData, ind::Int...) = extents.(Ref(dataset), ind)
cell_extents(dataset::AbstractVTKStructuredData) = ntuple(i -> extents(dataset, i) .- 1, dim(dataset))
cell_extents(dataset::AbstractVTKStructuredData, ind::Int) = cell_extents(dataset)[ind]
cell_extents(dataset::AbstractVTKStructuredData, ind::Int...) = cell_extents.(Ref(dataset), ind)

dim(dataset::AbstractVTKUnstructuredData) = size(dataset.point_coords, 1)
dim(dataset::AbstractVTKStructuredData) = length(extents(dataset))
function dim(dataset::AbstractVTKMultiblockData)
    return maximum(i -> dim(dataset[i]), 1:length(dataset))
end

num_of_points(dataset::AbstractVTKUnstructuredData) = size(dataset.point_coords, 2)
num_of_points(dataset::AbstractVTKStructuredData) = prod(extents(dataset))
num_of_points(dataset::AbstractVTKMultiblockData) = sum(num_of_points(i) for i in dataset)
num_of_points(dataset::AbstractTimeSeriesVTKData, ind::Int=1) = num_of_points(dataset[ind]) 

num_of_cells(dataset::AbstractVTKUnstructuredData) = length(dataset.cell_types)
num_of_cells(dataset::AbstractVTKStructuredData) = prod(extents(dataset) .- 1)
num_of_cells(dataset::AbstractVTKMultiblockData) = sum(num_of_cells(i) for i in dataset)
num_of_cells(dataset::AbstractTimeSeriesVTKData, ind::Int=1) = num_of_cells(dataset[ind])

num_of_point_vars(dataset::AbstractVTKSimpleData) = length(dataset.point_data)
num_of_cell_vars(dataset::AbstractVTKSimpleData) = length(dataset.cell_data)

function cell_type(dataset::Union{VTKUnstructuredData, VTKPolyData}, ind::Int)
    return dataset.cell_types[ind]
end
function cell_type(dataset::Union{VTKStructuredData, VTKRectilinearData}, ind::Int)
    if ind < 1 || ind > num_of_cells(dataset)
        throw("Out of bounds.")
    else
        return dim(dataset) == 2 ? 9 : 12
    end
end
function cell_type(dataset::VTKUniformRectilinearData, ind::Int)
    if ind < 1 || ind > num_of_cells(dataset)
        throw("Out of bounds.")
    else
        return dim(dataset) == 2 ? 8 : 11
    end
end
function cell_type(dataset::AbstractVTKStructuredData, ind::Tuple{Vararg{Int}})
    return cell_type(dataset, (LinearIndices(cell_extents(dataset)))[ind...])
end

function cell_connectivity(dataset::AbstractVTKUnstructuredData, cell_ind::Int)
    return dataset.cell_connectivity[cell_ind]
end
function cell_connectivity(dataset::T, cell_ind::Tuple{Vararg{Int}}) where {T<:AbstractVTKStructuredData}
    pextents = extents(dataset)
    return cell_connectivity(T, pextents, cell_ind)
end

function cell_connectivity(dataset::T, cell_ind) where {T<:AbstractVTKStructuredData}
    return cell_connectivity(T, extents(dataset), cell_ind)
end

function cell_connectivity(T::DataType, pextents::NTuple{N,Int}, cell_ind::NTuple{N,Int}) where {N}
    corner_point_ind = (LinearIndices(pextents))[cell_ind...]
    if N == 2 
        #1 <= cell_ind <= num_of_cells(dataset) || throw("Out of bounds.")
        if T <: VTKUniformRectilinearData
            return [corner_point_ind, corner_point_ind + 1, corner_point_ind + pextents[1], 
                    corner_point_ind + pextents[1] + 1]
        else
            return [corner_point_ind, corner_point_ind + 1, corner_point_ind + pextents[1] + 1, 
                    corner_point_ind + pextents[1]]
        end
    elseif N == 3
        #1 <= cell_ind <= num_of_cells(dataset) || throw("Out of bounds.")
        if T <: VTKUniformRectilinearData
            return [corner_point_ind, corner_point_ind + 1, corner_point_ind + pextents[1], 
                    corner_point_ind + pextents[1] + 1, 
                    corner_point_ind + pextents[1] * pextents[2], 
                    corner_point_ind + pextents[1] * pextents[2] + 1, 
                    corner_point_ind + pextents[1] * pextents[2] + pextents[1], 
                    corner_point_ind + pextents[1] * pextents[2] + pextents[1] + 1]
        else
            return [corner_point_ind, corner_point_ind + 1, corner_point_ind + pextents[1] + 1,          corner_point_ind + pextents[1], 
                    corner_point_ind + pextents[1] * pextents[2], 
                    corner_point_ind + pextents[1] * pextents[2] + 1, 
                    corner_point_ind + pextents[1] * pextents[2] + pextents[1] + 1, corner_point_ind + pextents[1] * pextents[2] + pextents[1]]
        end
    end

    throw("Invalid mesh dimension.")
end

function has_var(dataset::AbstractStaticVTKData, var_name::String)
    if haskey(dataset.point_data, var_name)
        return true, "Point"
    elseif haskey(dataset.cell_data, var_name)
        return true, "Cell"
    else
        return false, ""
    end
end

function var_dim(dataset::AbstractVTKUnstructuredData, var_name::String, var_type::String="")
    if var_type == "Point" && haskey(dataset.point_data, var_name)
        if length(size(dataset.point_data[var_name])) == 1
            return 1
        else
            return length(dataset.point_data[var_name][:,1])
        end
    elseif var_type == "Cell" && haskey(dataset.cell_data, var_name)
        if length(size(dataset.cell_data[var_name])) == 1
            return 1
        else
            return length(dataset.cell_data[var_name][:,1])
        end
    elseif var_type == ""
        try 
            return var_dim(dataset, var_name, "Point")
        catch
            return var_dim(dataset, var_name, "Cell")
        end
    else
        throw("Variable $var_name doesn't exist, please use has_var to check if the variable exists before using this function.")
    end
end

function var_dim(dataset::AbstractVTKStructuredData, var_name::String, var_type::String="")
    if var_type == "Point" && haskey(dataset.point_data, var_name)
        if length(size(dataset.point_data[var_name])) == length(extents(dataset))
            return 1
        else
            return size(dataset.point_data[var_name], 1)
        end
    elseif var_type == "Cell" && haskey(dataset.cell_data, var_name)
        if length(size(dataset.cell_data[var_name])) == length(extents(dataset))
            return 1
        else
            return size(dataset.cell_data[var_name], 1)
        end
    elseif var_type == ""
        try 
            return var_dim(dataset, var_name, "Point")
        catch
            return var_dim(dataset, var_name, "Cell")
        end
    else
        throw("Variable $var_name doesn't exist, please use has_var to check if the variable exists before using this function.")
    end
end

function vtk_cell_type_name(dataset::AbstractVTKStructuredData, cell_id)
    return VTK_CELL_TYPE[cell_type(dataset, cell_id)]
end
function is_homogeneous(dataset::AbstractVTKUnstructuredData)
    begin
        out = true
        for i in 2:num_of_cells(dataset)
            (out = cell_type(dataset, 1) == cell_type(dataset, i)) || break
        end
        out
    end || return false
    return true
end

is_homogeneous(dataset::AbstractVTKStructuredData) = true

function is_homogeneous(dataset::AbstractVTKMultiblockData)
    return all(i -> typeof(dataset[1]) == typeof(dataset[i]), 2:length(dataset))
end
function is_homogeneous(dataset::AbstractTimeSeriesVTKData)
    length(dataset) == 0 || 
    is_homogeneous(dataset[1]) && begin
        out = true
        i = dataset[1]
        T = typeof(i)
        for j in dataset
            (out = isa(j, T)) || break
            if T <: AbstractVTKMultiblockData || T <: AbstractVTKUnstructuredData 
                (out = same_ordered_geometry_shape(i,j)) || break
            else
                (out = same_geometry_shape(i,j)) || break
            end
        end
        out
    end || return false
    return true
end

function get_cell_ids(dataset::AbstractVTKUnstructuredData, cell_types::Vector{Int})
    cell_ids = Int[]
    for (i, c) in enumerate(dataset.cell_types)
        if in(c, cell_types)
            push!(cell_ids, i)
        end
    end
    return cell_ids
end

function get_lowest_index(cell_connectivity)
    return minimum(i -> minimum(cell_connectivity[i]), 1:length(cell_connectivity))
end
function get_highest_index(cell_connectivity)
    return maximum(i -> maximum(cell_connectivity[i]), 1:length(cell_connectivity))
end
function is_valid_cell(cell_connectivity, cell_type)
    return VTK_CELL_TYPE[cell_type].nodes == -1 || 
            length(_cell_connectivity) == VTK_CELL_TYPE[cell_type].nodes
end
num_of_blocks(dataset::AbstractVTKMultiblockData) = length(dataset.blocks)

timespan(dataset::VTKTimeSeriesData) = dataset.timemarkers[end] - dataset.timemarkers[1]

num_of_timesteps(dataset::VTKTimeSeriesData) = length(dataset.timemarkers)

triangular(dataset::AbstractVTKUnstructuredData) = all(i -> i == 5, dataset.cell_types)

function bb(dataset::T) where {T<:AbstractVTKData}
    if T <: AbstractVTKUnstructuredData
        return bb_unstruct(dataset)
    elseif T <: VTKStructuredData
        return bb_struct(dataset)
    elseif T <: VTKRectilinearData
        return bb_rect(dataset)
    elseif T <: VTKUniformRectilinearData
        return bb_image(dataset)
    elseif T <: Union{AbstractVTKMultiblockData, AbstractTimeSeriesVTKData}
        return bb_multiblock_time(dataset)
    end
    throw("Unsupported type.")
end

function bb_unstruct(dataset::AbstractVTKUnstructuredData)
    _dim = dim(dataset)
    min_x = max_x = dataset.point_coords[1,1]
    min_y = max_y = dataset.point_coords[2,1]
    if _dim == 2
        for i in 2:num_of_points(dataset)
            if dataset.point_coords[1,i] > max_x
                max_x = dataset.point_coords[1,i]
            elseif dataset.point_coords[1,i] < min_x
                min_x = dataset.point_coords[1,i]                
            end
        
            if dataset.point_coords[2,i] > max_y
                max_y = dataset.point_coords[2,i]
            elseif dataset.point_coords[2,i] < min_y
                min_y = dataset.point_coords[2,i]                
            end
        end
        return (min_x, max_x, min_y, max_y)
    elseif _dim == 3
        min_z = max_z = dataset.point_coords[3,1]
        for i in 2:num_of_points(dataset)
            if dataset.point_coords[1,i] > max_x
                max_x = dataset.point_coords[1,i]
            elseif dataset.point_coords[1,i] < min_x
                min_x = dataset.point_coords[1,i]                
            end
        
            if dataset.point_coords[2,i] > max_y
                max_y = dataset.point_coords[2,i]
            elseif dataset.point_coords[2,i] < min_y
                min_y = dataset.point_coords[2,i]                
            end

            if dataset.point_coords[3,i] > max_z
                max_z = dataset.point_coords[3,i]
            elseif dataset.point_coords[3,i] < min_z
                min_z = dataset.point_coords[3,i]                
            end
        end
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    else
        throw("Invalid dimension.")
    end
end

@resumable function surface_cell_inds(dataset::AbstractVTKStructuredData)
    cextents = cell_extents(dataset)
    if dim(dataset) == 2
        i = 1
        for j in 1:cextents[2]
            @yield (i,j)
        end
        j = 1
        for i in 2:cextents[1]
            @yield (i,j)
        end

        i = cextents[1]
        for j in 2:cextents[2]
            @yield (i,j)
        end
        j = cextents[2]
        for i in 2:cextents[1]-1
            @yield (i,j)
        end
    else
        i = 1
        for j in 1:cextents[2], k in 1:cextents[3]
            @yield (i,j,k)
        end
        
        j = 1
        for i in 2:cextents[1], k in 1:cextents[3]
            @yield (i,j,k)
        end

        k = 1
        for i in 2:cextents[1], j in 2:cextents[2]
            @yield (i,j,k)
        end

        i = cextents[1]
        for j in 2:cextents[2], k in 2:cextents[3]
            @yield (i,j,k)
        end
        
        j = cextents[2]
        for i in 2:cextents[1]-1, k in 2:cextents[3]
            @yield (i,j,k)
        end

        k = cextents[3]
        for i in 2:cextents[1]-1, j in 2:cextents[3]-1
            @yield (i,j,k)
        end
    end
end

function bb_struct(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _dim = dim(dataset)

    min_x = max_x = dataset.point_coords[1,1]
    min_y = max_y = dataset.point_coords[2,1]
    if _dim == 2
        for (i,j) in surface_cell_inds(dataset)
            if dataset.point_coords[1,i] > max_x
                max_x = dataset.point_coords[1,i]
            elseif dataset.point_coords[1,i] < min_x
                min_x = dataset.point_coords[1,i]                
            end
        
            if dataset.point_coords[2,i] > max_y
                max_y = dataset.point_coords[2,i]
            elseif dataset.point_coords[2,i] < min_y
                min_y = dataset.point_coords[2,i]                
            end
        end
        return (min_x, max_x, min_y, max_y)
    elseif _dim == 3
        min_z = max_z = dataset.point_coords[3,1]
        for (i,j,k) in surface_cell_inds(dataset)
            if dataset.point_coords[1,i] > max_x
                max_x = dataset.point_coords[1,i]
            elseif dataset.point_coords[1,i] < min_x
                min_x = dataset.point_coords[1,i]                
            end
        
            if dataset.point_coords[2,i] > max_y
                max_y = dataset.point_coords[2,i]
            elseif dataset.point_coords[2,i] < min_y
                min_y = dataset.point_coords[2,i]                
            end

            if dataset.point_coords[3,i] > max_z
                max_z = dataset.point_coords[3,i]
            elseif dataset.point_coords[3,i] < min_z
                min_z = dataset.point_coords[3,i]                
            end
        end
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    else
        throw("Invalid dimension.")
    end
end

function bb_rect(dataset::VTKRectilinearData)
    min_x = dataset.point_coords[1][1]
    max_x = dataset.point_coords[1][end]
    min_y = dataset.point_coords[2][1]
    max_y = dataset.point_coords[2][end]

    if dim(dataset) == 2
        return (min_x, max_x, min_y, max_y)
    elseif dim(dataset) == 3
        min_z = dataset.point_coords[3][1]
        max_z = dataset.point_coords[3][end]
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    else
        throw("Invalid dimension.")
    end
end

function bb_image(dataset::VTKUniformRectilinearData)
    min_x = dataset.origin[1]
    max_x = dataset.origin[1] + (dataset.extents[1]-1)*dataset.spacing[1]
    min_x = dataset.origin[2]
    max_x = dataset.origin[2] + (dataset.extents[2]-1)*dataset.spacing[2]

    if dim(dataset) == 2
        return (min_x, max_x, min_y, max_y)
    elseif dim(dataset) == 3
        min_x = dataset.origin[3]
        max_x = dataset.origin[3] + (dataset.extents[3]-1)*dataset.spacing[3]
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    else
        throw("Invalid dimension.")
    end
end

function bb_multiblock_time(dataset::Union{AbstractVTKMultiblockData, AbstractTimeSeriesVTKData})
    _bbs = [bb(block) for block in dataset]
    min_x = minimum(i -> _bbs[i][1], 1:length(_bbs))
    max_x = maximum(i -> _bbs[i][2], 1:length(_bbs))
    min_y = minimum(i -> _bbs[i][3], 1:length(_bbs))
    max_y = maximum(i -> _bbs[i][4], 1:length(_bbs))

    _first = true
    min_z = max_z = 0
    for _bb in _bbs
        if length(_bb) == 6
            if _first
                min_z = _bb[5]
                max_z = _bb[6]
                _first = false
            end
            if _bb[5] < min_z
                min_z = _bb[5]
            end
            if _bb[5] > max_z
                max_z = _bb[6]
            end
        end
    end

    if min_z == max_z == 0
        return (min_x, max_x, min_y, max_y)
    else
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    end
end

function pseudo_center(dataset::T) where {T<:AbstractVTKData}
    _bb = bb(dataset)
    return ntuple(i -> (_bb[2i-1] + _bb[2i])/2, Val(length(_bb)÷2))
end
