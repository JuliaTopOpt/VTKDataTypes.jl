import Base.convert

convert{T<:AbstractVTKData}(::Type{T}, dataset::T) = dataset

convert(::Type{VTKMultiblockData}, dataset::AbstractVTKSimpleData) = VTKMultiblockData(AbstractStaticVTKData[dataset])

function convert(::Type{VTKRectilinearData}, dataset::VTKUniformRectilinearData)
    _extents = extents(dataset)
    _dim = dim(dataset)
    point_coords = [[dataset.origin[j] + (i-1)*dataset.spacing[j] for i in 1:_extents[j]] for j in 1:_dim]
    return VTKRectilinearData(point_coords, dataset.point_data, dataset.cell_data)
end

convert(::Type{VTKStructuredData}, dataset::VTKUniformRectilinearData) = VTKStructuredData(VTKRectilinearData(dataset))

function convert(::Type{VTKStructuredData}, dataset::VTKRectilinearData)
    if dim(dataset) == 2
        xcoord = dataset.point_coords[1]
        ycoord = dataset.point_coords[2]

        point_coords = zeros(dim(dataset), extents(dataset)...)
        for j in 1:extents(dataset, 2), i in 1:extents(dataset, 1)
            point_coords[:, i, j] = [xcoord[i], ycoord[j]]
        end
    elseif dim(dataset) == 3
        xcoord = dataset.point_coords[1]
        ycoord = dataset.point_coords[2]
        zcoord = dataset.point_coords[3]

        point_coords = zeros(dim(dataset), extents(dataset)...)
        for k in 1:extents(dataset, 3), j in 1:extents(dataset, 2), i in 1:extents(dataset, 1)
            point_coords[:, i, j, k] = [xcoord[i], ycoord[j], zcoord[k]]
        end
    else
        throw("Invalid dataset dimensions.")
    end
    return VTKStructuredData(point_coords, dataset.point_data, dataset.cell_data)
end

function convert(::Type{VTKUnstructuredData}, dataset::VTKPolyData)
    return VTKUnstructuredData(dataset.point_coords, dataset.cell_types, 
        dataset.cell_connectivity, dataset.point_data, dataset.cell_data)
end

function convert{T<:AbstractVTKStructuredData}(::Type{VTKUnstructuredData}, dataset::T)
    point_coords = zeros(dim(dataset), num_of_points(dataset))
    if T <: VTKStructuredData
        point_coords = unstructured_point_coords_from_structured(dataset)
    elseif T <: VTKRectilinearData
        point_coords = unstructured_point_coords_from_rectilinear(dataset)
    elseif T <: VTKUniformRectilinearData
        point_coords = unstructured_point_coords_from_uniform_rectilinear(dataset)
    end

    point_data = unstructured_point_data_from_structured(dataset)
    cell_types, _cell_connectivity, cell_data = unstructured_cell_info_from_structured(dataset)
    return VTKUnstructuredData(point_coords, cell_types, _cell_connectivity, point_data, cell_data)
end

function convert(::Type{VTKUnstructuredData}, data_blocks::VTKMultiblockData)
    combined_block = VTKUnstructuredData(deepcopy(data_blocks[1]))
    for block in data_blocks[2:end]
        append!(combined_block, VTKUnstructuredData(block))
    end
    return combined_block
end

convert(::Type{VTKPolyData}, dataset::AbstractStaticVTKData) = 
    decompose(VTKUnstructuredData(dataset), "Faces", true)

function unstructured_point_data_from_structured(dataset::AbstractVTKStructuredData)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)
    _extents = extents(dataset)

    point_data = Dict{String, Array{Float64}}()
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (num_of_points(dataset),))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, num_of_points(dataset)))
        end
    end
    point_data
end

function unstructured_cell_info_from_structured{T<:AbstractVTKStructuredData}(dataset::T)
    _num_of_cells = num_of_cells(dataset)
    _dim = dim(dataset)
    pextents = extents(dataset)
    cextents = cell_extents(dataset)

    local_cell_connectivity(cell_ind) = cell_connectivity(T, pextents, cell_ind)
    if _dim == 2
        if T <: VTKUniformRectilinearData
            cell_types = fill(8, _num_of_cells) #VTK_PIXEL
        else
            cell_types = fill(9, _num_of_cells) #VTK_QUAD
        end
        _cell_connectivity = reshape(map(local_cell_connectivity, ((i,j) for j in 1:cextents[2], i in 1:cextents[1])), (_num_of_cells,))
    elseif _dim == 3
        if T <: VTKUniformRectilinearData
            cell_types = fill(11, _num_of_cells) #VTK_VOXEL
        else
            cell_types = fill(12, _num_of_cells) #VTK_HEXAHEDRON
        end
        _cell_connectivity = reshape(map(local_cell_connectivity, ((i,j,k) for k in 1:cextents[3], j in 1:cextents[2], i in 1:cextents[1])), (_num_of_cells,))
    end

    cell_data = Dict{String, Array{Float64}}()
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            cell_data[m] = reshape(dataset.cell_data[m], (num_of_cells(dataset),))
        else
            cell_data[m] = reshape(dataset.cell_data[m], (_var_dim, num_of_cells(dataset)))
        end
    end
    cell_types, _cell_connectivity, cell_data
end

function unstructured_point_coords_from_structured(dataset::VTKStructuredData)
    _extents = extents(dataset)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)

    return reshape(dataset.point_coords, (_dim, _num_of_points))
end

function unstructured_point_coords_from_rectilinear(dataset::VTKRectilinearData)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)

    _point_coords = zeros(_dim, _num_of_points)
    p = 1
    if _dim == 2
        for y in dataset.point_coords[2], x in dataset.point_coords[1]
            _point_coords[1,p] = x
            _point_coords[2,p] = y
            p += 1
        end
    elseif _dim == 3
        for z in dataset.point_coords[3], y in dataset.point_coords[2], x in dataset.point_coords[1]
            _point_coords[1,p] = x
            _point_coords[2,p] = y
            _point_coords[3,p] = z
            p += 1
        end
    else
        throw("Invalid dimension.")
    end
    return _point_coords
end

function unstructured_point_coords_from_uniform_rectilinear(dataset::VTKUniformRectilinearData)
    _extents = extents(dataset)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)
    _spacing = dataset.spacing
    _origin = dataset.origin

    _point_coords = zeros(_dim, _num_of_points)
    p = 1
    if _dim == 2
        x, y = _origin
        for j in 1:_extents[2]
            x = _origin[1]
            for i in 1:_extents[1]
                _point_coords[1,p] = x
                _point_coords[2,p] = y
                p += 1
                x += _spacing[1]
            end
            y += _spacing[2]
        end
    elseif _dim == 3
        x, y, z = _origin
        for k in 1:_extents[3]
            y = _origin[2]
            for j in 1:_extents[2]
                x = _origin[1]
                for i in 1:_extents[1]
                    _point_coords[1,p] = x
                    _point_coords[2,p] = y
                    _point_coords[3,p] = z
                    p += 1
                    x += _spacing[1]
                end
                y += _spacing[2]
            end
            z += _spacing[3]
        end
    else
        throw("Invalid dimension.")
    end

    return _point_coords
end
