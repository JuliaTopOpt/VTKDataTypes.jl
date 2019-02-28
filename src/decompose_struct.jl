function decompose(_dataset::AbstractVTKStructuredData, target::String = "Faces", decompose_cell_data = false)
    dataset = VTKStructuredData(_dataset)
    _dim = dim(dataset)

    if decompose_cell_data
        if target == "Faces"
            if _dim == 2
                return decompose_to_faces_2d_with_cell_data(dataset)
            elseif _dim == 3
                return decompose_to_faces_3d_with_cell_data(dataset)
            else
                throw("Unsupported dimension.")
            end
        elseif target == "Lines"
            if _dim == 2
                return decompose_to_lines_2d_with_cell_data(dataset)
            elseif _dim == 3
                return decompose_to_lines_3d_with_cell_data(dataset)
            else
                throw("Unsupported dimension.")
            end
        else
            throw("Unsupported target.")
        end
    else
        if target == "Faces"
            if _dim == 2
                return decompose_to_faces_2d_no_cell_data(dataset)
            elseif _dim == 3
                return decompose_to_faces_3d_no_cell_data(dataset)
            else
                throw("Unsupported dimension.")
            end
        elseif target == "Lines"
            if _dim == 2
                return decompose_to_lines_2d_no_cell_data(dataset)
            elseif _dim == 3
                return decompose_to_lines_3d_no_cell_data(dataset)
            else
                throw("Unsupported dimension.")
            end
        else
            throw("Unsupported target.")
        end
    end
end

function decompose_to_faces_2d_with_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)
    _num_of_cells = num_of_cells(dataset)
    _dim = dim(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    faces = Vector{Int}[]
    _cell_types = fill(9, _num_of_cells)
    for cind in Iterators.product(1:cextents[1], 1:cextents[2])
        push!(faces, cell_connectivity(dataset, cind))
    end
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            _cell_data[m] = reshape(dataset.cell_data[m], (_num_of_cells,))
        else
            _cell_data[m] = reshape(dataset.cell_data[m], (_var_dim, _num_of_cells))
        end
    end
    return VTKPolyData(point_coords, _cell_types, faces, point_data, _cell_data)
end

function decompose_to_faces_3d_with_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    faces = Vector{Int}[]
    for m in keys(dataset.cell_data)
        _cell_data[m] = Float64[]
    end

    for cind in Iterators.product(1:cextents[1], 1:cextents[2], 1:cextents[3])
        hexa_cc = cell_connectivity(dataset, cind)
        push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[3], hexa_cc[4]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[3] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1:2]...,cind[3]-1]) / 2
                else
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[3] != 1
                    _cd = (dataset.cell_data[m][:, cind...] +
                        dataset.cell_data[m][:, cind[1:2]...,cind[3]-1]) / 2
                else
                    _cd = dataset.cell_data[m][:, cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[6], hexa_cc[5]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[2] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1],cind[2]-1,cind[3]]) / 2
                else
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[3] != 1
                    _cd = (dataset.cell_data[m][:, cind...] +
                        dataset.cell_data[m][:, cind[1],cind[2]-1,cind[3]]) / 2
                else
                    _cd = dataset.cell_data[m][:, cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        push!(faces, [hexa_cc[4], hexa_cc[1], hexa_cc[5], hexa_cc[8]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[1] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1]-1,cind[2:3]...]) / 2
                else
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[3] != 1
                    _cd = (dataset.cell_data[m][:, cind...] +
                        dataset.cell_data[m][:, cind[1]-1,cind[2:3]...]) / 2
                else
                    _cd = dataset.cell_data[m][:, cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        if cind[1] == cextents[1]
            push!(faces, [hexa_cc[2], hexa_cc[3], hexa_cc[7], hexa_cc[6]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    _cd = dataset.cell_data[m][cind...]
                    push!(_cell_data[m], _cd)
                else
                    _cd = dataset.cell_data[m][:, cind...]
                    append!(_cell_data[m], _cd)
                end
            end
        end
        if cind[2] == cextents[2]
            push!(faces, [hexa_cc[4], hexa_cc[3], hexa_cc[7], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    _cd = dataset.cell_data[m][cind...]
                    push!(_cell_data[m], _cd)
                else
                    _cd = dataset.cell_data[m][:, cind...]
                    append!(_cell_data[m], _cd)
                end
            end
        end
        if cind[3] == cextents[3]
            push!(faces, [hexa_cc[5], hexa_cc[6], hexa_cc[7], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    _cd = dataset.cell_data[m][cind...]
                    push!(_cell_data[m], _cd)
                else
                    _cd = dataset.cell_data[m][:, cind...]
                    append!(_cell_data[m], _cd)
                end
            end
        end
    end
    _cell_types = fill(9, length(faces))

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            _cell_data[m] = reshape(_cell_data[m], (length(faces),))
        else
            _cell_data[m] = reshape(_cell_data[m], (_var_dim, length(faces)))
        end
    end

    return VTKPolyData(point_coords, _cell_types, faces, point_data, _cell_data)
end

function decompose_to_lines_2d_with_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    lines = Vector{Int}[]
    for m in keys(dataset.cell_data)
        _cell_data[m] = Float64[]
    end
    for cind in Iterators.product(1:cextents[1], 1:cextents[2])
        quad_cc = cell_connectivity(dataset, cind)

        push!(lines, [quad_cc[1], quad_cc[2]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[2] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1], cind[2]-1]) / 2
                else
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[2] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1],cind[2]-1]) / 2
                else
                    _cd = dataset.cell_data[m][:,cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        push!(lines, [quad_cc[1], quad_cc[4]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[1] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1]-1, cind[2]]) / 2
                else
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[1] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1]-1,cind[2]]) / 2
                else
                    _cd = dataset.cell_data[m][:,cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        if cind[1] == cextents[1]
            push!(lines, [quad_cc[2], quad_cc[3]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    _cd = dataset.cell_data[m][cind...]
                    push!(_cell_data[m], _cd)
                else
                    _cd = dataset.cell_data[m][:,cind...]
                    append!(_cell_data[m], _cd)
                end
            end
        end
        if cind[2] == cextents[2]
            push!(lines, [quad_cc[3], quad_cc[4]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    _cd = dataset.cell_data[m][cind...]
                    push!(_cell_data[m], _cd)
                else
                    _cd = dataset.cell_data[m][:,cind...]
                    append!(_cell_data[m], _cd)
                end
            end
        end
        _cell_types = fill(3, length(lines))
    end

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            _cell_data[m] = reshape(dataset.cell_data[m], (length(lines),))
        else
            _cell_data[m] = reshape(dataset.cell_data[m], (_var_dim, length(lines)))
        end
    end

    return VTKPolyData(point_coords, _cell_types, lines, point_data, _cell_data)
end

function decompose_to_lines_3d_with_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    lines = Vector{Int}[]
    for m in keys(dataset.cell_data)
        _cell_data[m] = Float64[]
    end

    for cind in Iterators.product(1:cextents[1], 1:cextents[2], 1:cextents[3])
        hexa_cc = cell_connectivity(dataset, cind)
        push!(lines, [hexa_cc[1], hexa_cc[2]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[2] != 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1:2]...,cind[3]-1] + 
                        dataset.cell_data[m][cind[1],cind[2]-1,cind[3]] + 
                        dataset.cell_data[m][cind[1],cind[2]-1,cind[3]-1]) / 4                            
                elseif cind[2] == 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1:2]...,cind[3]-1]) / 2
                elseif cind[2] != 1 && cind[3] == 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1], cind[2]-1,cind[3]]) / 2                            
                elseif cind[2] == cind[3] == 1
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[2] != 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1:2]...,cind[3]-1] + 
                        dataset.cell_data[m][:,cind[1],cind[2]-1,cind[3]] + 
                        dataset.cell_data[m][:,cind[1],cind[2]-1,cind[3]-1]) / 4                            
                elseif cind[2] == 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1:2]...,cind[3]-1]) / 2
                elseif cind[2] != 1 && cind[3] == 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1], cind[2]-1,cind[3]]) / 2                            
                elseif cind[2] == cind[3] == 1
                    _cd = dataset.cell_data[m][:,cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        push!(lines, [hexa_cc[1], hexa_cc[4]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[1] != 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1:2]...,cind[3]-1] + 
                        dataset.cell_data[m][cind[1]-1,cind[2],cind[3]] + 
                        dataset.cell_data[m][cind[1]-1,cind[2],cind[3]-1]) / 4                            
                elseif cind[1] == 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1:2]...,cind[3]-1]) / 2
                elseif cind[1] != 1 && cind[3] == 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1]-1, cind[2],cind[3]]) / 2                            
                elseif cind[1] == cind[3] == 1
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[1] != 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1:2]...,cind[3]-1] + 
                        dataset.cell_data[m][:,cind[1]-1,cind[2],cind[3]] + 
                        dataset.cell_data[m][:,cind[1]-1,cind[2],cind[3]-1]) / 4                            
                elseif cind[1] == 1 && cind[3] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1:2]...,cind[3]-1]) / 2
                elseif cind[1] != 1 && cind[3] == 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1]-1, cind[2],cind[3]]) / 2                            
                elseif cind[1] == cind[3] == 1
                    _cd = dataset.cell_data[m][:,cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        push!(lines, [hexa_cc[1], hexa_cc[5]])
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                if cind[1] != 1 && cind[2] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1]-1,cind[2:3]...] + 
                        dataset.cell_data[m][cind[1],cind[2]-1,cind[3]] + 
                        dataset.cell_data[m][cind[1]-1,cind[2]-1,cind[3]]) / 4                            
                elseif cind[1] == 1 && cind[2] != 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1],cind[2]-1,cind[3]]) / 2
                elseif cind[1] != 1 && cind[2] == 1
                    _cd = (dataset.cell_data[m][cind...] +
                        dataset.cell_data[m][cind[1]-1,cind[2],cind[3]]) / 2
                elseif cind[1] == cind[2] == 1
                    _cd = dataset.cell_data[m][cind...]
                end
                push!(_cell_data[m], _cd)
            else
                if cind[1] != 1 && cind[2] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1]-1,cind[2:3]...] + 
                        dataset.cell_data[m][:,cind[1],cind[2]-1,cind[3]] + 
                        dataset.cell_data[m][:,cind[1]-1,cind[2]-1,cind[3]]) / 4                            
                elseif cind[1] == 1 && cind[2] != 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1],cind[2]-1,cind[3]]) / 2
                elseif cind[1] != 1 && cind[2] == 1
                    _cd = (dataset.cell_data[m][:,cind...] +
                        dataset.cell_data[m][:,cind[1]-1,cind[2],cind[3]]) / 2
                elseif cind[1] == cind[2] == 1
                    _cd = dataset.cell_data[m][:,cind...]
                end
                append!(_cell_data[m], _cd)
            end
        end

        if cind[1] == cextents[1]
            push!(lines, [hexa_cc[2], hexa_cc[3]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[3] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1:2]..., cind[3]-1]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[3] != 1
                        _cd = (dataset.cell_data[m][:,cind...] +
                            dataset.cell_data[m][:,cind[1:2]...,cind[3]-1]) / 2
                    else
                        _cd = dataset.cell_data[m][:,cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[3], hexa_cc[7]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[2] != cextents[2]
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1], cind[2]+1, cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[2] != cextents[2]
                        _cd = (dataset.cell_data[m][:,cind...] +
                            dataset.cell_data[m][:,cind[1],cind[2]+1,cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][:,cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[6], hexa_cc[7]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[3] != cextents[3]
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1:2]..., cind[3]+1]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[3] != cextents[3]
                        _cd = (dataset.cell_data[m][:,cind...] +
                            dataset.cell_data[m][:,cind[1:2]...,cind[3]+1]) / 2
                    else
                        _cd = dataset.cell_data[m][:,cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[2], hexa_cc[6]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[2] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1], cind[2]-1, cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[2] != 1
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1], cind[2]-1, cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end
        end
        if cind[2] == cextents[2]
            push!(lines, [hexa_cc[3], hexa_cc[4]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[3] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1:2]..., cind[3]-1]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[3] != 1
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1:2]..., cind[3]-1]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[7], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[3] != cextents[3]
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1:2]..., cind[3]+1]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[3] != cextents[3]
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1:2]..., cind[3]+1]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[4], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[1] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1]-1, cind[2:3]...]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[1] != 1
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1]-1, cind[2:3]...]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end
        end
        if cind[3] == cextents[3]
            push!(lines, [hexa_cc[5], hexa_cc[8]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[1] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1]-1, cind[2:3]...]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[1] != 1
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1]-1, cind[2:3]...]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end

            push!(lines, [hexa_cc[5], hexa_cc[6]])
            for m in keys(dataset.cell_data)
                _var_dim = var_dim(dataset, m, "Cell")
                if _var_dim == 1
                    if cind[2] != 1
                        _cd = (dataset.cell_data[m][cind...] +
                            dataset.cell_data[m][cind[1], cind[2]-1, cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][cind...]
                    end
                    push!(_cell_data[m], _cd)
                else
                    if cind[2] != 1
                        _cd = (dataset.cell_data[m][:, cind...] +
                            dataset.cell_data[m][:, cind[1], cind[2]-1, cind[3]]) / 2
                    else
                        _cd = dataset.cell_data[m][:, cind...]
                    end
                    append!(_cell_data[m], _cd)
                end
            end
        end
    end
    _cell_types = fill(3, length(lines))

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            _cell_data[m] = reshape(dataset.cell_data[m], (length(lines),))
        else
            _cell_data[m] = reshape(dataset.cell_data[m], (_var_dim, length(lines)))
        end
    end

    return VTKPolyData(point_coords, _cell_types, lines, point_data, _cell_data)
end

function decompose_to_faces_2d_no_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)
    _num_of_cells = num_of_cells(dataset)
    _dim = dim(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    faces = Vector{Int}[]
    _cell_types = fill(9, _num_of_cells)
    for cind in Iterators.product(1:cextents[1], 1:cextents[2])
        push!(faces, cell_connectivity(dataset, cind))
    end
    return VTKPolyData(point_coords, _cell_types, faces, point_data, _cell_data)
end

function decompose_to_faces_3d_no_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)
    _dim = dim(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    faces = Vector{Int}[]

    for cind in Iterators.product(1:cextents[1], 1:cextents[2], 1:cextents[3])
        hexa_cc = cell_connectivity(dataset, cind)
        push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[3], hexa_cc[4]])
        push!(faces, [hexa_cc[1], hexa_cc[2], hexa_cc[6], hexa_cc[5]])
        push!(faces, [hexa_cc[4], hexa_cc[1], hexa_cc[5], hexa_cc[8]])
        if cind[1] == cextents[1]
            push!(faces, [hexa_cc[2], hexa_cc[3], hexa_cc[7], hexa_cc[6]])
        end
        if cind[2] == cextents[2]
            push!(faces, [hexa_cc[4], hexa_cc[3], hexa_cc[7], hexa_cc[8]])
        end
        if cind[3] == cextents[3]
            push!(faces, [hexa_cc[5], hexa_cc[6], hexa_cc[7], hexa_cc[8]])
        end
    end
    _cell_types = fill(9, length(faces))

    return VTKPolyData(point_coords, _cell_types, faces, point_data, _cell_data)
end

function decompose_to_lines_2d_no_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    lines = Vector{Int}[]

    for cind in Iterators.product(1:cextents[1], 1:cextents[2])
        quad_cc = cell_connectivity(dataset, cind)

        push!(lines, [quad_cc[1], quad_cc[2]])
        push!(lines, [quad_cc[1], quad_cc[4]])
        if cind[1] == cextents[1]
            push!(lines, [quad_cc[2], quad_cc[3]])
        end
        if cind[2] == cextents[2]
            push!(lines, [quad_cc[3], quad_cc[4]])
        end
        _cell_types = fill(3, length(lines))
    end

    return VTKPolyData(point_coords, _cell_types, lines, point_data, _cell_data)
end

function decompose_to_lines_3d_no_cell_data(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    _num_of_points = num_of_points(dataset)

    point_coords = reshape(dataset.point_coords, (_dim, _num_of_points))
    point_data = empty(dataset.point_data)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            point_data[m] = reshape(dataset.point_data[m], (_num_of_points,))
        else
            point_data[m] = reshape(dataset.point_data[m], (_var_dim, _num_of_points))
        end
    end

    _cell_data = empty(dataset.cell_data)
    lines = Vector{Int}[]

    for cind in Iterators.product(1:cextents[1], 1:cextents[2], 1:cextents[3])
        hexa_cc = cell_connectivity(dataset, cind)
        push!(lines, [hexa_cc[1], hexa_cc[2]])
        push!(lines, [hexa_cc[1], hexa_cc[4]])
        push!(lines, [hexa_cc[1], hexa_cc[5]])
        if cind[1] == cextents[1]
            push!(lines, [hexa_cc[2], hexa_cc[3]])
            push!(lines, [hexa_cc[3], hexa_cc[7]])
            push!(lines, [hexa_cc[6], hexa_cc[7]])
            push!(lines, [hexa_cc[2], hexa_cc[6]])
        end
        if cind[2] == cextents[2]
            push!(lines, [hexa_cc[3], hexa_cc[4]])
            push!(lines, [hexa_cc[7], hexa_cc[8]])
            push!(lines, [hexa_cc[4], hexa_cc[8]])
        end
        if cind[3] == cextents[3]
            push!(lines, [hexa_cc[5], hexa_cc[8]])
            push!(lines, [hexa_cc[5], hexa_cc[6]])
        end
    end
    _cell_types = fill(3, length(lines))

    return VTKPolyData(point_coords, _cell_types, lines, point_data, _cell_data)
end
