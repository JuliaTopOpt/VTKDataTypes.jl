
function extract_surface(_dataset::AbstractVTKStructuredData)
    if dim(_dataset) == 2
        return VTKPolyData(_dataset)
    end

    dataset = VTKStructuredData(_dataset)
    cextents = cell_extents(dataset)
    faces = Vector{Int}[]

    _cell_data = empty(_dataset.cell_data)
    for m in keys(dataset.cell_data)
        _cell_data[m] = Float64[]
    end

    for cind in surface_cell_inds(dataset)
        hexa_cc = cell_connectivity(dataset, cind)
        if cind[1] == 1
            push!(faces, [hexa_cc[1], hexa_cc[4], hexa_cc[8], hexa_cc[5]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end
        if cind[1] == cextents[1]
            push!(faces, [hexa_cc[2], hexa_cc[3], hexa_cc[7], hexa_cc[6]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end
        
        if cind[2] == 1
            push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[6], hexa_cc[5]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end
        if cind[2] == cextents[2]
            push!(faces, [hexa_cc[3], hexa_cc[4], hexa_cc[8], hexa_cc[7]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end

        if cind[3] == 1
            push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[3], hexa_cc[4]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end
        if cind[3] == cextents[3]
            push!(faces, [hexa_cc[5], hexa_cc[6], hexa_cc[7], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                add_cell_data!(dataset, _cell_data, m, cind)
            end
        end
    end

    _cell_types = fill(9, length(faces))

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim > 1
            _cell_data[m] = reshape(_cell_data[m], (_var_dim, length(faces)))
        end
    end

    _point_coords = Float64[]
    point_inds = Int[]
    _point_data = empty(_dataset.point_data)
    for m in keys(dataset.point_data)
        _point_data[m] = Float64[]
    end
    pextents = extents(dataset)

    k = 0
    for _face in faces
        for (i,p) in enumerate(_face)
            _i = findfirst(point_inds, p)
            if _i > 0
                _face[i] = _i
            else
                k += 1
                append!(_point_coords, dataset.point_coords[:,ind2sub(pextents, p)...])
                push!(point_inds, p)
                for m in keys(dataset.point_data)
                    _var_dim = var_dim(dataset, m, "Point")
                    if _var_dim == 1
                        push!(_point_data[m], dataset.point_data[m][ind2sub(pextents, p)...])
                    else
                        append!(_point_data[m], dataset.point_data[m][:, ind2sub(pextents, p)...])
                    end
                end
                _face[i] = k
            end
        end
    end

    _point_coords = reshape(_point_coords, (3, k))
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            _point_data[m] = reshape(_point_data[m], (k,))
        else
            _point_data[m] = reshape(_point_data[m], (_var_dim, k))
        end
    end

    return VTKPolyData(_point_coords, _cell_types, faces, _point_data, _cell_data)
    #=
    _point_coords = reshape(dataset.point_coords, (3, num_of_points(dataset)))
    _point_data = empty(_dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            _point_data[m] = reshape(dataset.point_data[m], (num_of_points(dataset),))            
        else
            _point_data[m] = reshape(dataset.point_data[m], (_var_dim, num_of_points(dataset)))
        end
    end
    return VTKPolyData(_point_coords, _cell_types, faces, _point_data, _cell_data)
    =#
end

function add_cell_data!(dataset, _cell_data, m, cind)
    _var_dim = var_dim(dataset, m, "Cell")
    if _var_dim == 1
        _cd = dataset.cell_data[m][cind...]
        push!(_cell_data[m], _cd)
    else
        _cd = dataset.cell_data[m][:, cind...]
        append!(_cell_data[m], _cd)
    end
    return
end

function extract_surface(dataset::T) where {T<:AbstractVTKUnstructuredData}
    if T <: VTKPolyData
        return dataset
    end
    if dim(dataset) == 2
        return VTKPolyData(dataset)
    end
    
    #Sorted inds and cell type are keys
    #Cell connectivity, cell count, and cell data are the values

    cell_register = Dict{Tuple{Vector{Int}, Int}, Tuple{Vector{Int}, _Counter, typeof(dataset.cell_data)}}()
    for i in 1:length(dataset.cell_connectivity)
        if dataset.cell_types[i] ∉ POLY_CELLS
            _cells, _types = decompose_cell(dataset.cell_connectivity[i], dataset.cell_types[i], target="Faces")
        else
            _cells, _types = [dataset.cell_connectivity[i]], [dataset.cell_types[i]]
        end
        for j in 1:length(_cells)
            _key = (sort(_cells[j]), _types[j])
            if haskey(cell_register, _key)
                cell_register[_key][2].a += 1
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

    _cell_data = empty(dataset.cell_data)
    _cell_connectivity = Vector{Int}[]
    _cell_types = Int[]

    for m in keys(dataset.cell_data)
        _cell_data[m] = Float64[]
    end

    for (i, kv) in enumerate(cell_register)
        k, v = kv
        if v[2].a == 1
            push!(_cell_connectivity, v[1])
            push!(_cell_types, k[2])
            for m in keys(dataset.cell_data)
                append!(_cell_data[m], v[3][m])
            end
        end
    end

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim > 1
            _cell_data[m] = reshape(_cell_data[m], (_var_dim, length(_cell_connectivity)))
        end
    end

    _point_coords = Float64[]
    point_inds = Int[]
    _point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _point_data[m] = Float64[]
    end

    k = 0
    for _cc in _cell_connectivity
        for (i,p) in enumerate(_cc)
            _i = findfirst(point_inds, p)
            if _i > 0
                _cc[i] = _i
            else
                k += 1
                append!(_point_coords, dataset.point_coords[:,p])
                push!(point_inds, p)
                for m in keys(dataset.point_data)
                    _var_dim = var_dim(dataset, m, "Point")
                    if _var_dim == 1
                        push!(_point_data[m], dataset.point_data[m][p])
                    else
                        append!(_point_data[m], dataset.point_data[m][:,p])
                    end
                end
                _cc[i] = k
            end
        end
    end

    _point_coords = reshape(_point_coords, (3, k))
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim > 1
            _point_data[m] = reshape(_point_data[m], (_var_dim, k))
        end
    end

    return VTKPolyData(_point_coords, _cell_types, _cell_connectivity, _point_data, _cell_data)
end
