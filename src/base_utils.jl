import Base: promote_rule, size, length, getindex, iterate, ==, append!, +, *, /

promote_rule(::Type{<:VTKUniformRectilinearData}, ::Type{<:VTKRectilinearData}) = VTKRectilinearData
promote_rule(::Type{<:VTKUniformRectilinearData}, ::Type{<:VTKStructuredData}) = VTKStructuredData
promote_rule(::Type{<:VTKRectilinearData}, ::Type{<:VTKStructuredData}) = VTKStructuredData
promote_rule(::Type{<:VTKUniformRectilinearData}, ::Type{<:VTKUnstructuredData}) = VTKUnstructuredData
promote_rule(::Type{<:VTKRectilinearData}, ::Type{<:VTKUnstructuredData}) = VTKUnstructuredData
promote_rule(::Type{<:VTKStructuredData}, ::Type{<:VTKUnstructuredData}) = VTKUnstructuredData
promote_rule(::Type{<:VTKMultiblockData}, ::Type{T}) where {T <: AbstractVTKSimpleData} = VTKMultiblockData
promote_rule(::Type{T}, ::Type{T}) where {T<:AbstractVTKData} = T

size(dataset::AbstractTimeSeriesVTKData) = (num_of_timesteps(dataset),)
size(dataset::AbstractVTKMultiblockData) = (num_of_blocks(dataset),)
length(dataset::Union{AbstractVTKMultiblockData, AbstractTimeSeriesVTKData}) = prod(size(dataset))

getindex(dataset::AbstractTimeSeriesVTKData, ind::Integer) = dataset.data[ind]
getindex(dataset::AbstractVTKMultiblockData, ind::Integer) = dataset.blocks[ind]
getindex(dataset::AbstractTimeSeriesVTKData, inds::OrdinalRange{Int,Int}) = dataset.data[inds]
function getindex(dataset::AbstractTimeSeriesVTKData, t::AbstractFloat)
    loc = searchsorted(dataset.timemarkers, t)
    if loc.start == 1
        return dataset[1]
    elseif loc.stop == length(dataset)
        return dataset[end]
    else
        return ((dataset.timemarkers[loc.start] - t) * dataset[loc.stop] + 
                (t - dataset.timemarkers[loc.stop]) * dataset[loc.start]) / 
                (dataset.timemarkers[loc.start] - dataset.timemarkers[loc.stop])
    end
end
getindex(dataset::AbstractVTKMultiblockData, inds::OrdinalRange{Int,Int}) = dataset.blocks[inds]

function iterate(dataset::AbstractTimeSeriesVTKData, i=1)
    i > length(dataset) && return nothing
    return dataset.data[i], i+1
end
function iterate(dataset::AbstractVTKMultiblockData, i=1)
    i > length(dataset) && return nothing
    return dataset.blocks[i], i+1
end

function +(a::T, b::T) where {T<:AbstractStaticVTKData}
    if !(T <: AbstractVTKMultiblockData) || !(same_ordered_geometry_shape(a,b))
        same_geometry_shape(a,b) || throw("Cannot add two datasets with different geometry structure.")
    end
    if !(T <: AbstractVTKMultiblockData)
        same_data_shape(a,b) || throw("Cannot add two datasets with different data variables.")
    end
    if T <: AbstractVTKMultiblockData
        return VTKMultiblockData(ntuple(i -> a[i] + b[i], Val(length(a))))
    elseif T <: AbstractVTKUnstructuredData || T <: AbstractVTKStructuredData
        point_coords = a.point_coords + b.point_coords
    elseif T <: VTKRectilinearData
        point_coords = [a.point_coords[i] + b.point_coords[i] for i in 1:dim(a)]
    elseif T <: VTKUniformRectilinearData
        origin = a.origin + b.origin
        spacing = a.spacing + b.spacing
    end

    point_data = empty(a.point_data)
    for m in keys(a.point_data)
        point_data[m] = a.point_data[m] + b.point_data[m]
    end
    cell_data = empty(a.cell_data)
    for m in keys(a.cell_data)
        cell_data[m] = a.cell_data[m] + b.cell_data[m]
    end

    if T <: AbstractVTKUnstructuredData
        return T(point_coords, a.cell_types, a.cell_connectivity, point_data, cell_data)
    elseif T <: AbstractVTKStructuredData || T <: VTKRectilinearData
        return T(point_coords, point_data, cell_data)
    elseif T <: VTKUniformRectilinearData
        return T(origin, spacing, a.extents, point_data, cell_data)
    end
end

*(b::Real, a::AbstractStaticVTKData) = *(a,b)
function *(a::T, b::Real) where {T <: AbstractStaticVTKData}
    if T <: AbstractVTKMultiblockData
        return VTKMultiblockData(ntuple(i -> *(a[i], b), Val(length(a))))
    elseif T <: AbstractVTKUnstructuredData || T <: AbstractVTKStructuredData
        point_coords = a.point_coords .* b
    elseif T <: VTKRectilinearData
        point_coords = [a.point_coords[i] .* b for i in 1:dim(a)]
    elseif T <: VTKUniformRectilinearData
        origin = a.origin .* b
        spacing = a.spacing .* b
    end

    point_data = empty(a.point_data)
    for m in keys(a.point_data)
        point_data[m] = a.point_data[m] .* b
    end
    cell_data = empty(a.cell_data)
    for m in keys(a.cell_data)
        cell_data[m] = a.cell_data[m] .* b
    end

    if T <: AbstractVTKUnstructuredData
        return T(point_coords, a.cell_types, a.cell_connectivity, point_data, cell_data)
    elseif T <: AbstractVTKStructuredData || T <: VTKRectilinearData
        return T(point_coords, point_data, cell_data)
    elseif T <: VTKUniformRectilinearData
        return T(origin, spacing, a.extents, point_data, cell_data)
    end
end

function /(a::T, b::Real) where {T <: AbstractStaticVTKData}
    b != 0 || throw("Cannot divide a dataset by zero.")

    if T <: AbstractVTKMultiblockData
        return VTKMultiblockData(ntuple(i -> a[i] / b, Val(length(a))))
    elseif T <: AbstractVTKUnstructuredData || T <: AbstractVTKStructuredData
        point_coords = a.point_coords ./ b
    elseif T <: VTKRectilinearData
        point_coords = [a.point_coords[i] ./ b for i in 1:dim(a)]
    elseif T <: VTKUniformRectilinearData
        origin = a.origin ./ b
        spacing = a.spacing ./ b
    end

    point_data = empty(a.point_data)
    for m in keys(a.point_data)
        point_data[m] = a.point_data[m] ./ b
    end
    cell_data = empty(a.cell_data)
    for m in keys(a.cell_data)
        cell_data[m] = a.cell_data[m] ./ b
    end

    if T <: AbstractVTKUnstructuredData
        return T(point_coords, a.cell_types, a.cell_connectivity, point_data, cell_data)
    elseif T <: AbstractVTKStructuredData || T <: VTKRectilinearData
        return T(point_coords, point_data, cell_data)
    elseif T <: VTKUniformRectilinearData
        return T(origin, spacing, a.extents, point_data, cell_data)
    end
end

function ==(a::T,b::T) where {T<:AbstractVTKStructuredData}
    return same_geometry(a,b) && a.point_data == b.point_data && a.cell_data == b.cell_data
end
function ==(a::T, b::T) where {T<:VTKUniformRectilinearData}
    return  dim(a) == dim(b) && extents(a) == extents(b) && 
            num_of_points(a) == num_of_points(b) && num_of_cells(a) == num_of_cells(b) && 
            a.origin == b.origin && a.spacing == b.spacing && a.point_data == b.point_data && 
            a.cell_data == b.cell_data
end
function ==(a::T, b::T) where {T<:VTKMultiblockData}
    return length(a) == length(b) && a ⊆ b && b ⊆ a
end
function ==(a::AbstractVTKUnstructuredData, b::AbstractVTKUnstructuredData)
    dim(a) != dim(b) && return false
    num_of_points(a) != num_of_points(b) && return false
    num_of_cells(a) != num_of_cells(b) && return false
    
    _num_of_points = num_of_points(a)
    #Describing mapping between ith cell in b and jth cell in a
    b_image_point_inds = zeros(Int, _num_of_points)
    for i in 1:_num_of_points
        for j in 1:_num_of_points
            @views if b.point_coords[:,i] == a.point_coords[:,j]
                b_image_point_inds[i] = j
                break
            end
        end
        if b_image_point_inds[i] == 0
            return false
        end
    end

    #Describing a cell of b in terms of points of a
    function map_cell_connectivity(b_cell_connectivity)
        b_cell_connectivity_mapped = zeros(Int, length(b_cell_connectivity))
        for i in 1:length(b_cell_connectivity)
            b_cell_connectivity_mapped[i] = b_image_point_inds[b_cell_connectivity[i]]
        end
        b_cell_connectivity_mapped
    end

    #Mapping the ith cell of b and the jth cell of a
    _num_of_cells = num_of_cells(a)
    b_cell_connectivity_mapped = map(map_cell_connectivity, b.cell_connectivity)
    b_image_cell_inds = zeros(Int, _num_of_cells)
    for i in 1:_num_of_cells
        for j in 1:_num_of_cells
            if similar_cells( b_cell_connectivity_mapped[i], 
                              cell_type(b,i), 
                              a.cell_connectivity[j], 
                              cell_type(a,j)
                            )
                b_image_cell_inds[i] = j
                break
            end
        end
        if b_image_cell_inds[i] == 0
            return false
        end
    end    

    a_keys = keys(a.point_data)
    b_keys = keys(b.point_data)
    length(a_keys) != length(b_keys) && return false
    for m in a_keys
        if !in(m, b_keys)
            return false
        end
    end
    a_keys = keys(a.cell_data)
    b_keys = keys(b.cell_data)
    length(a_keys) != length(b_keys) && return false
    for m in a_keys
        if !in(m, b_keys)
            return false
        end
    end

    for v in keys(a.point_data)
        _var_dim = var_dim(a, v, "Point")
        if _var_dim == 1
            for i in 1:_num_of_points
                b.point_data[v][i] != a.point_data[v][b_image_point_inds[i]] && return false
            end        
        else
            @views for i in 1:_num_of_points
                b.point_data[v][:,i] != a.point_data[v][:,b_image_point_inds[i]] && return false
            end
        end            
    end

    for v in keys(a.cell_data)
        _var_dim = var_dim(a, v, "Cell")
        if _var_dim == 1
            for i in 1:_num_of_cells
                b.cell_data[v][i] == a.cell_data[v][b_image_cell_inds[i]] || return false
            end
        else
            @views for i in 1:_num_of_cells
                b.cell_data[v][:,i] == a.cell_data[v][:,b_image_cell_inds[i]] || return false
            end
        end
    end

    return true
end

#=
function ==(a::AbstractVTKUnstructuredData, b::AbstractVTKUnstructuredData)
    dim(a) != dim(b) && return false
    num_of_points(a) != num_of_points(b) && return false
    num_of_cells(a) != num_of_cells(b) && return false
    
    _num_of_points = num_of_points(a)
    _num_of_cells = num_of_cells(a)
    #Describing mapping between ith cell in b and jth cell in a
    
    a_keys = keys(a.point_data)
    b_keys = keys(b.point_data)
    length(a_keys) != length(b_keys) && return false
    for m in a_keys
        if !in(m, b_keys)
            return false
        end
    end
    pkeys = a_keys

    for i in 1:_num_of_points
        check = false
        for j in 1:_num_of_points
            if a.point_coords[:,i] == b.point_coords[:,j]
                for m in pkeys
                    if var_dim(a, m, "Point") == 1
                        a.point_data[m][i] == b.point_data[m][j] || return false
                    else
                        a.point_data[m][:,i] == b.point_data[m][:,j] || return false
                    end
                end
                check = true
                break
            end
        end
        if !check
            return false
        end
    end

    a_keys = keys(a.cell_data)
    b_keys = keys(b.cell_data)
    length(a_keys) != length(b_keys) && return false
    for m in a_keys
        if !in(m, b_keys)
            return false
        end
    end
    ckeys = a_keys
    for i in 1:_num_of_cells
        check = false
        for j in 1:_num_of_cells
            a.cell_types[i] == b.cell_types[j] || continue
            points_matched = 0
            for a_p in a.cell_connectivity[i]
                check2 = false
                for b_p in b.cell_connectivity[j]
                    if a.point_coords[:,a_p] == b.point_coords[:,b_p]
                        check2 = true
                        points_matched += 1
                        break
                    end
                end
                if !check2
                    break
                end
            end
            if points_matched != length(a.cell_connectivity[i])
                continue
            end

            for m in ckeys
                if var_dim(a, m, "Cell") == 1
                    a.cell_data[m][i] == b.cell_data[m][j] || return false
                else
                    a.cell_data[m][:,i] == b.cell_data[m][:,j] || return false
                end
            end
            check = true
            break
        end
        if !check
            return false
        end
    end
    return true
end
=#

function ==(a::T, b::T) where {T<:AbstractTimeSeriesVTKData}
    length(a) == length(b) || return false
    a.timemarkers == b.timemarkers || return false
    a.data == b.data || return false
    return true
end

==(a::AbstractVTKData, b::AbstractVTKData) = false

function append!(datasets::AbstractVTKUnstructuredData...)
    dataset1 = datasets[1]
    for i in 2:length(datasets)
        point_offset = num_of_points(dataset1)
        cell_offset = num_of_cells(dataset1)
        dataset1.point_coords = [dataset1.point_coords datasets[i].point_coords]
        dataset1.cell_types = [dataset1.cell_types; datasets[i].cell_types]
        dataset1.cell_connectivity = [dataset1.cell_connectivity; datasets[i].cell_connectivity]
        for j in (cell_offset + 1):length(dataset1.cell_connectivity)
            dataset1.cell_connectivity[j] .+= point_offset
        end
        
        for m in keys(dataset1.point_data)
            _var_dim = var_dim(dataset1, m, "Point")
            if in(m, keys(datasets[i].point_data))
                if _var_dim == 1
                    dataset1.point_data[m] = [dataset1.point_data[m]; datasets[i].point_data[m]]
                else
                    dataset1.point_data[m] = [dataset1.point_data[m] datasets[i].point_data[m]]
                end
            end
        end

        for m in keys(dataset1.cell_data)
            _var_dim = var_dim(dataset1, m, "Cell")
            if in(m, keys(datasets[i].cell_data))
                if _var_dim == 1
                    dataset1.cell_data[m] = [dataset1.cell_data[m]; datasets[i].cell_data[m]]
                else
                    dataset1.cell_data[m] = [dataset1.cell_data[m] datasets[i].cell_data[m]]
                end
            end
        end
    end

    return dataset1
end
