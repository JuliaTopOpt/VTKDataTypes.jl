
function triangulate_with_cell_data{S<:Real}(dataset::AbstractVTKUnstructuredData{S})
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity, cell count, and cell data are the values

    cell_register = Dict{Vector{Int}, Tuple{Vector{Int}, _Counter, Dict{String, Array{Float64}}}}()
    for i in 1:length(dataset.cell_connectivity)
        _cells = triangulate_cell(dataset.cell_connectivity[i], dataset.cell_types[i])
        for j in 1:length(_cells)
            _key = sort(_cells[j])
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
                cell_register[_key] = (_cells[j], _Counter(1), Dict{String, Array{Float64}}())
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
    cell_data = Dict{String, Array{Float64}}()
    _cell_connectivity = Vector{Int}[]

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
        push!(_cell_connectivity, v[1])
    end

    _cell_types = [5 for i in 1:ncells]

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function triangulate_no_cell_data{S<:Real}(dataset::AbstractVTKUnstructuredData{S})
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity is the value
    cell_register = Dict{Vector{Int}, Vector{Int}}()
    for i in 1:length(dataset.cell_connectivity)
        _cells = triangulate_cell(dataset.cell_connectivity[i], dataset.cell_types[i])
        for j in 1:length(_cells)
            _key = sort(_cells[j])
            if !haskey(cell_register, _key)
                cell_register[_key] = _cells[j]
            end
        end
    end

    ncells = length(cell_register)
    cell_data = Dict{String, Array{Float64}}()
    _cell_connectivity = Vector{Int}[]

    for (i, kv) in enumerate(cell_register)
        k, v = kv
        push!(_cell_connectivity, v)
    end
    _cell_types = [5 for i in 1:ncells]

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function triangulate(dataset::AbstractVTKUnstructuredData, triangulate_cell_data=false)
    if triangulate_cell_data
        return triangulate_with_cell_data(dataset)
    else
        return triangulate_no_cell_data(dataset)
    end
end

