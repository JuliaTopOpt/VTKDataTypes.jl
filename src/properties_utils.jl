function similar_cells(cell_connectivity1::Vector{Int}, cell_type1::Int, 
    cell_connectivity2::Vector{Int}, cell_type2::Int)
    for k in cell_connectivity2
        in(k, cell_connectivity1) || return false
    end
    cell_type1 == cell_type2 || return false
    return true
end

coherent{T<:AbstractVTKMultiblockData}(a::T, b::T) = same_ordered_geometry_shape(a,b) && same_data_shape(a,b)
coherent{T<:AbstractVTKSimpleData}(a::T, b::T) = same_geometry_shape(a, b) && same_data_shape(a, b)
function coherent{T<:AbstractTimeSeriesVTKData}(dataset::T)
    length(dataset) == 0 && return true 
    block1 = dataset[1]
    for block in dataset[2:end]
        coherent(block1, block) || return false
    end
    return true
end
function coherent{T<:AbstractTimeSeriesVTKData}(a::T, b::T)
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
coherent{T<:AbstractVTKData, S<:AbstractVTKData}(a::T, b::S) = false

extents(dataset::VTKStructuredData) = size(dataset.point_coords)[2:end]
extents(dataset::VTKRectilinearData) = ([length(dataset.point_coords[i]) for i in 1:length(dataset.point_coords)]...)
extents(dataset::VTKUniformRectilinearData) = (dataset.extents...)
extents(dataset::AbstractVTKStructuredData, ind::Int) = extents(dataset)[ind]
extents(dataset::AbstractVTKStructuredData, ind::Int...) = extents(dataset)[[ind...]]
cell_extents(dataset::AbstractVTKStructuredData) = (([extents(dataset)...] .- 1)...)
cell_extents(dataset::AbstractVTKStructuredData, ind::Int) = cell_extents(dataset)[ind]
cell_extents(dataset::AbstractVTKStructuredData, ind::Int...) = cell_extents(dataset)[[ind...]]

dim(dataset::AbstractVTKUnstructuredData) = size(dataset.point_coords, 1)
dim(dataset::AbstractVTKStructuredData) = length(extents(dataset))
function dim(dataset::AbstractVTKMultiblockData)
    _dims = Int[]
    for b in dataset
        push!(_dims, dim(b))
    end
    return max(_dims...)
end

num_of_points(dataset::AbstractVTKUnstructuredData) = size(dataset.point_coords, 2)
num_of_points(dataset::AbstractVTKStructuredData) = reduce(*, extents(dataset))
num_of_points(dataset::AbstractVTKMultiblockData) = sum([num_of_points(i) for i in dataset])
num_of_points(dataset::AbstractTimeSeriesVTKData, ind::Int=1) = num_of_points(dataset[ind]) 

num_of_cells(dataset::AbstractVTKUnstructuredData) = length(dataset.cell_types)
num_of_cells(dataset::AbstractVTKStructuredData) = reduce(*, [extents(dataset)...].-1)
num_of_cells(dataset::AbstractVTKMultiblockData) = sum([num_of_cells(i) for i in dataset])
num_of_cells(dataset::AbstractTimeSeriesVTKData, ind::Int=1) = num_of_cells(dataset[ind]) 

num_of_point_vars(dataset::AbstractVTKSimpleData) = length(dataset.point_data)
num_of_cell_vars(dataset::AbstractVTKSimpleData) = length(dataset.cell_data)

cell_type(dataset::Union{VTKUnstructuredData, VTKPolyData}, ind::Int) = dataset.cell_types[ind]
cell_type(dataset::Union{VTKStructuredData, VTKRectilinearData}, ind::Int) = ind < 1 || ind > num_of_cells(dataset) ? throw("Out of bounds.") : dim(dataset) == 2 ? 9 : 12
cell_type(dataset::VTKUniformRectilinearData, ind::Int) = ind < 1 || ind > num_of_cells(dataset) ? throw("Out of bounds.") : dim(dataset) == 2 ? 8 : 11

cell_connectivity(dataset::AbstractVTKUnstructuredData, cell_ind::Int) = dataset.cell_connectivity[cell_ind]
function cell_connectivity{T<:AbstractVTKStructuredData, N<:Integer}(dataset::T, _cell_ind::NTuple{N, Int})
    pextents = extents(dataset)
    cell_connectivity(T, pextents, _cell_ind)
end

function cell_connectivity{T<:AbstractVTKStructuredData}(dataset::T, cell_ind)
    return cell_connectivity(T, extents(dataset), cell_ind)
end

function cell_connectivity{N}(T::DataType, pextents::NTuple{N,Int}, _cell_ind::NTuple{N,Int})
    corner_point_ind = sub2ind(pextents, _cell_ind...)
    if N == 2 
        #1 <= cell_ind <= num_of_cells(dataset) || throw("Out of bounds.")
        if T <: VTKUniformRectilinearData
            return [corner_point_ind, corner_point_ind+1, corner_point_ind+pextents[1], corner_point_ind+pextents[1]+1]
        else
            return [corner_point_ind, corner_point_ind+1, corner_point_ind+pextents[1]+1, corner_point_ind+pextents[1]]
        end
    elseif N == 3
        #1 <= cell_ind <= num_of_cells(dataset) || throw("Out of bounds.")
        if T <: VTKUniformRectilinearData
            return [corner_point_ind, corner_point_ind+1, corner_point_ind+pextents[1], corner_point_ind+pextents[1]+1, corner_point_ind+pextents[1]*pextents[2], corner_point_ind+pextents[1]*pextents[2]+1, corner_point_ind+pextents[1]*pextents[2]+pextents[1], corner_point_ind+pextents[1]*pextents[2]+pextents[1]+1]
        else
            return [corner_point_ind, corner_point_ind+1, corner_point_ind+pextents[1]+1, corner_point_ind+pextents[1], corner_point_ind+pextents[1]*pextents[2], corner_point_ind+pextents[1]*pextents[2]+1, corner_point_ind+pextents[1]*pextents[2]+pextents[1]+1, corner_point_ind+pextents[1]*pextents[2]+pextents[1]]
        end
    end

    throw("Invalid mesh dimension.")
end

function has_var(dataset::AbstractStaticVTKData, var_name::String)
    if haskey(dataset.point_data, var_name)
        true, "Point"
    elseif haskey(dataset.cell_data, var_name)
        true, "Cell"
    else
        false, ""
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

vtk_cell_type_name(dataset::AbstractVTKStructuredData, cell_id) = VTK_CELL_TYPE[cell_type(dataset, cell_id)]

function is_homogeneous(dataset::AbstractVTKUnstructuredData)
    begin
        _out = true
        for i in 2:num_of_cells(dataset)
            (_out = cell_type(dataset, 1) == cell_type(dataset, i)) || break
        end
        _out
    end || return false
    return true
end

is_homogeneous(dataset::AbstractVTKStructuredData) = true

is_homogeneous(dataset::AbstractVTKMultiblockData) = all([typeof(dataset[1]) == typeof(dataset[i]) for i in 2:length(dataset)])

function is_homogeneous(dataset::AbstractTimeSeriesVTKData)
    length(dataset) == 0 || 
    is_homogeneous(dataset[1]) && begin
        _out = true
        i = dataset[1]
        T = typeof(i)
        for j in dataset
            (_out = isa(j, T)) || break
            if T <: AbstractVTKMultiblockData || T <: AbstractVTKUnstructuredData 
                (_out = same_ordered_geometry_shape(i,j)) || break
            else
                (_out = same_geometry_shape(i,j)) || break
            end
        end
        _out
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
    cell_ids
end

get_lowest_index(_cell_connectivity::Vector{Vector{Int}}) = min([min(i...) for i in _cell_connectivity]...)

get_highest_index(_cell_connectivity::Vector{Vector{Int}}) = max([max(i...) for i in _cell_connectivity]...)

is_valid_cell(_cell_connectivity::Vector{Int}, cell_type) = VTK_CELL_TYPE[cell_type].nodes == -1 || 
    length(_cell_connectivity) == VTK_CELL_TYPE[cell_type].nodes

num_of_blocks(dataset::AbstractVTKMultiblockData) = length(dataset.blocks)

timespan(dataset::VTKTimeSeriesData) = dataset.timemarkers[end] - dataset.timemarkers[1]

num_of_timesteps(dataset::VTKTimeSeriesData) = length(dataset.timemarkers)

triangular(dataset::AbstractVTKUnstructuredData) = all(dataset.cell_types .== 5)

function bb{T<:AbstractVTKData}(dataset::T)
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
        return [min_x, max_x, min_y, max_y]
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
        return [min_x, max_x, min_y, max_y, min_z, max_z]
    else
        throw("Invalid dimension.")
    end
end

@pygen function surface_cell_inds(dataset::AbstractVTKStructuredData)
    cextents = cell_extents(dataset)
    if dim(dataset) == 2
        i = 1
        for j in 1:cextents[2]
            yield((i,j))
        end
        j = 1
        for i in 2:cextents[1]
            yield((i,j))
        end

        i = cextents[1]
        for j in 2:cextents[2]
            yield((i,j))
        end
        j = cextents[2]
        for i in 2:cextents[1]-1
            yield((i,j))
        end
    else
        i = 1
        for j in 1:cextents[2], k in 1:cextents[3]
            yield((i,j,k))
        end
        
        j = 1
        for i in 2:cextents[1], k in 1:cextents[3]
            yield((i,j,k))
        end

        k = 1
        for i in 2:cextents[1], j in 2:cextents[2]
            yield((i,j,k))
        end

        i = cextents[1]
        for j in 2:cextents[2], k in 2:cextents[3]
            yield((i,j,k))
        end
        
        j = cextents[2]
        for i in 2:cextents[1]-1, k in 2:cextents[3]
            yield((i,j,k))
        end

        k = cextents[3]
        for i in 2:cextents[1]-1, j in 2:cextents[3]-1
            yield((i,j,k))
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
        return [min_x, max_x, min_y, max_y]
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
        return [min_x, max_x, min_y, max_y, min_z, max_z]
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
        return [min_x, max_x, min_y, max_y]
    elseif dim(dataset) == 3
        min_z = dataset.point_coords[3][1]
        max_z = dataset.point_coords[3][end]
        return [min_x, max_x, min_y, max_y, min_z, max_z]
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
        return [min_x, max_x, min_y, max_y]
    elseif dim(dataset) == 3
        min_x = dataset.origin[3]
        max_x = dataset.origin[3] + (dataset.extents[3]-1)*dataset.spacing[3]
        return [min_x, max_x, min_y, max_y, min_z, max_z]
    else
        throw("Invalid dimension.")
    end
end

function bb_multiblock_time(dataset::Union{AbstractVTKMultiblockData, AbstractTimeSeriesVTKData})
    _bbs = [bb(block) for block in dataset]
    min_x = min([_bbs[i][1] for i in 1:length(_bbs)]...)
    max_x = max([_bbs[i][2] for i in 1:length(_bbs)]...)
    min_y = min([_bbs[i][3] for i in 1:length(_bbs)]...)
    max_y = max([_bbs[i][4] for i in 1:length(_bbs)]...)

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
        return [min_x, max_x, min_y, max_y]
    else
        return [min_x, max_x, min_y, max_y, min_z, max_z]
    end
end

function pseudo_center{T<:AbstractVTKData}(dataset::T)
    _bb = bb(dataset)
    return [(_bb[i] + _bb[i+1])/2 for i in 1:2:length(_bb)]
end
