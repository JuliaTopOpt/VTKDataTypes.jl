function triangulate_with_cell_data(dataset::AbstractVTKUnstructuredData, OutputT = NTuple{3,Int})
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity, cell count, and cell data are the values

    cell_register = Dict{NTuple{3, Int}, Tuple{NTuple{3, Int}, _Counter, typeof(dataset.cell_data)}}()
    for i in 1:length(dataset.cell_connectivity)
        _cells = triangulate_cell(dataset.cell_connectivity[i], dataset.cell_types[i])
        for j in 1:length(_cells)
            _key = sort(SVector{3,Int}(_cells[j])).data
            if haskey(cell_register, _key)
                cell_register[_key][2].a += 1
                for m in keys(dataset.cell_data)
                    _var_dim = var_dim(dataset, m, "Cell")
                    if _var_dim == 1
                        cell_register[_key][3][m] += [dataset.cell_data[m][i]]
                    else
                        cell_register[_key][3][m] += dataset.cell_data[m][:,i]
                    end
                end
            else
                cell_register[_key] = (_cells[j], _Counter(1), empty(dataset.cell_data))
                for m in keys(dataset.cell_data)
                    _var_dim = var_dim(dataset, m, "Cell")
                    if _var_dim == 1
                        cell_register[_key][3][m] = [dataset.cell_data[m][i]]
                    else
                        cell_register[_key][3][m] = dataset.cell_data[m][:,i]
                    end
                end
            end
        end
    end

    ncells = length(cell_register)
    cell_data = empty(dataset.cell_data)
    _cell_connectivity = OutputT[]

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            cell_data[m] = zeros(length(cell_register))
        else
            cell_data[m] = zeros(_var_dim, length(cell_register))
        end
    end

    for (i, kv) in enumerate(cell_register)
        k, v = kv
        for m in keys(dataset.cell_data)
            v[3][m] = v[3][m] ./ v[2].a
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                cell_data[m][i] = v[3][m][1]
            else
                cell_data[m][:,i] = v[3][m]
            end
        end
        push!(_cell_connectivity, OutputT(v[1]...))
    end

    _cell_types = [5 for i in 1:ncells]

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function triangulate_no_cell_data(dataset::AbstractVTKUnstructuredData, OutputT = NTuple{3, Int})
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity is the value
    cell_register = Dict{NTuple{3, Int}, NTuple{3, Int}}()
    for i in 1:length(dataset.cell_connectivity)
        _cells = triangulate_cell(dataset.cell_connectivity[i], dataset.cell_types[i])
        for j in 1:length(_cells)
            _key = sort(SVector{3,Int}(_cells[j])).data
            if !haskey(cell_register, _key)
                cell_register[_key] = _cells[j]
            end
        end
    end

    ncells = length(cell_register)
    cell_data = empty(dataset.cell_data)
    _cell_connectivity = OutputT[]

    for (i, kv) in enumerate(cell_register)
        k, v = kv
        push!(_cell_connectivity, OutputT(v...))
    end
    _cell_types = [5 for i in 1:ncells]

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function duplicate_vertices(dataset::AbstractVTKUnstructuredData)
    @assert all(isequal(5), dataset.cell_types)
    point_coords = zeros(size(dataset.point_coords, 1), 3*num_of_cells(dataset))
    cell_connectivity = deepcopy(dataset.cell_connectivity)
    point_data = empty(dataset.point_data)
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            point_data[m] = zeros(3*num_of_cells(dataset))
        else
            point_data[m] = zeros(_var_dim, 3*num_of_points(dataset))
        end
    end

    point_counter = 0
    for i in 1:num_of_cells(dataset)
        for j in dataset.cell_connectivity[i]
            point_counter += 1
            @views point_coords[:,point_counter] += dataset.point_coords[:,j]
            for m in keys(dataset.point_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    point_data[m][point_counter] += dataset.point_data[m][j]
                else
                    @views point_data[m][:,point_counter] += dataset.point_data[m][:,j]
                end
            end
        end
        cell_connectivity[i] = (point_counter-2, point_counter-1, point_counter)
    end

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            point_data[m] = point_data[m] ./ max.(point_counter, 1)
        else
            point_data[m] = point_data[m] ./ max.(point_counter', 1)
        end
    end
    cell_data = deepcopy(dataset.cell_data)
    cell_types = copy(dataset.cell_types)

    return VTKPolyData(point_coords, cell_types, cell_connectivity, point_data, cell_data)
end

function triangulate(dataset::AbstractVTKUnstructuredData, triangulate_cell_data=false, OutputT = NTuple{3,Int})
    if triangulate_cell_data
        return triangulate_with_cell_data(dataset, OutputT)
    else
        return triangulate_no_cell_data(dataset, OutputT)
    end
end
