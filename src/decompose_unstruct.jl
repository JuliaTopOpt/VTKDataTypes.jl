
function decompose(dataset::AbstractVTKUnstructuredData, target::String="Faces", decompose_cell_data=false)
    if decompose_cell_data
        return decompose_with_cell_data(dataset, target)
    else
        return decompose_no_cell_data(dataset, target)
    end
end

type _Counter
    a::Int
end

function decompose_with_cell_data{S<:Real}(dataset::AbstractVTKUnstructuredData{S}, target::String="Faces")
    if target == "Faces"
        filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])
    elseif target == "Lines"
        filter_cells!(dataset, POINT_CELLS)
    end

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity, cell count, and cell data are the values

    cell_register = Dict{Tuple{Vector{Int}, Int}, Tuple{Vector{Int}, _Counter, Dict{String, Array}}}()
    for i in 1:length(dataset.cell_connectivity)
        if target == "Points"
            _cells, _types = [[k] for k in dataset.cell_connectivity[i]], [1 for k in 1:length(dataset.cell_connectivity[i])]
        else
            _cells, _types = decompose_cell(dataset.cell_connectivity[i], dataset.cell_types[i], target=target)
        end        
        for j in 1:length(_cells)
            _key = (sort(_cells[j]), _types[j])
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
                cell_register[_key] = (_cells[j], _Counter(1), Dict{String, Array}())
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
    _cell_connectivity = [Int[] for i in 1:ncells]
    _cell_types = zeros(Int, ncells)

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
        _cell_connectivity[i] = v[1]
        _cell_types[i] = k[2]
    end

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function decompose_no_cell_data{S<:Real}(dataset::AbstractVTKUnstructuredData{S}, target::String="Faces")
    if target == "Faces"
        filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])
    elseif target == "Lines"
        filter_cells!(dataset, POINT_CELLS)
    end

    point_coords = dataset.point_coords
    point_data = dataset.point_data
    
    #Sorted inds and cell type are keys
    #Cell connectivity is the value
    cell_register = Dict{Tuple{Vector{Int}, Int}, Vector{Int}}()
    for i in 1:length(dataset.cell_connectivity)
        if target == "Points"
            _cells, _types =  [[k] for k in dataset.cell_connectivity[i]], [1 for k in 1:length(dataset.cell_connectivity[i])]
        else
            _cells, _types = decompose_cell(dataset.cell_connectivity[i], dataset.cell_types[i], target=target)
        end        
        for j in 1:length(_cells)
            _key = (sort(_cells[j]), _types[j])
            if !haskey(cell_register, _key)
                cell_register[_key] = _cells[j]
            end
        end
    end

    ncells = length(cell_register)
    cell_data = Dict{String, Array{Float64}}()
    _cell_connectivity = [Int[] for i in 1:ncells]
    _cell_types = Int[]

    for (i, kv) in enumerate(cell_register)
        k, v = kv
        _cell_connectivity[i] = v
        push!(_cell_types, k[2])
    end

    return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end
