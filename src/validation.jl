function is_valid_unstructured_mesh{_T<:AbstractVTKUnstructuredData}(dataset::_T, repeat_cells::Bool)
    T = repr(_T)
    _point_coords, _cell_types, _cell_connectivity, _point_data, _cell_data = dataset.point_coords, 
        dataset.cell_types, dataset.cell_connectivity, dataset.point_data, dataset.cell_data

    in(dim(dataset), [2,3]) || return false, "$T: Only 2 and 3 dimensional mesh data are allowed."

    length(_cell_types) == length(_cell_connectivity) || return false, "$T: The cell types vector and the cell connectivity vector have unequal lengths."
    length(_cell_connectivity) == 0 || (_min = get_lowest_index(_cell_connectivity)) >= 1 || return false, "$T: The lowest point index used in the cell connectivity vector is $_min < 1. Please consider using the add_point_id_offset!(cell_connectivity, $(1-_min)) function to make the lowest index 1."
    length(_cell_connectivity) == 0 || get_highest_index(_cell_connectivity) <= num_of_points(dataset) || return false, "$T: The highest point index used in the cell connectivity vector is more than the number of points passed."

    _keys = collect(keys(_point_data))
    length(_keys) == 0 || 
    num_of_points(dataset) == (var_dim(dataset, _keys[1], "Point") == 1 ? length(_point_data[_keys[1]]) : size(_point_data[_keys[1]], 2)) && 
    begin 
        _out = true
        _var_dim1 = var_dim(dataset, _keys[1], "Point")
        for m in _keys[2:end]
            _var_dimi = var_dim(dataset, m, "Point")
            (_out = (_var_dim1 == 1 ? length(_point_data[_keys[1]]) : size(_point_data[_keys[1]], 2)) == 
                (_var_dimi == 1 ? length(_point_data[m]) : size(_point_data[m], 2))) || break
        end
        _out
    end || return false, "$T: The number of points defined is not consistent with the number of points in at least one of the point data fields."

    _keys = collect(keys(_cell_data))
    length(_keys) == 0 || 
    num_of_cells(dataset) == (var_dim(dataset, _keys[1], "Cell") == 1 ? length(_cell_data[_keys[1]]) : size(_cell_data[_keys[1]], 2)) && 
    begin 
        _out = true
        _var_dim1 = var_dim(dataset, _keys[1], "Cell")
        for m in _keys[2:end]
            _var_dimi = var_dim(dataset, m, "Cell")
            (_out = (_var_dim1 == 1 ? length(_cell_data[_keys[1]]) : size(_cell_data[_keys[1]], 2)) == 
                (_var_dimi == 1 ? length(_cell_data[m]) : size(_cell_data[m], 2))) || break
        end
        _out
    end || return false, "$T: The number of cells defined is not consistent with the number of cells in at least one of the cell data fields."

    begin
        _out = true
        for i in 1:length(_cell_connectivity)
            (_out = VTK_CELL_TYPE[cell_type(dataset, i)].nodes == -1 || 
                length(_cell_connectivity[i]) == VTK_CELL_TYPE[cell_type(dataset, i)].nodes) || break
        end
        _out
    end || return false, "$T: The number of points in at least one cell does not match its corresponding cell type."

    if !repeat_cells
        begin
            _out = true
            for i in 1:length(_cell_connectivity)
                for j in i+1:length(_cell_connectivity)
                    (_out = !similar_cells(_cell_connectivity[i], cell_type(dataset, i), 
                        _cell_connectivity[j], cell_type(dataset, j))) || break
                end
            end
            _out
        end || return false, "$T: Repeat cells are not allowed."
    end
    
    return true, ""
end

function is_valid(dataset::VTKUnstructuredData; repeat_cells=false)
    return is_valid_unstructured_mesh(dataset, repeat_cells)
end

function is_valid(dataset::VTKPolyData; repeat_cells=false)
    T = "VTKPolyData"
    _cell_types, _cell_connectivity = dataset.cell_types, dataset.cell_connectivity
    
    length(_cell_types) == 0 || begin 
        _out = true
        poly_cells = [POINT_CELLS; LINE_CELLS; FACE_CELLS]
        for i in 1:length(_cell_connectivity)
            (_out = in(cell_type(dataset, i), poly_cells)) || break
        end
        _out
    end || return false, "$T: At least one of the cells is a volume cell. $T can only have point, line or face cells, consider VTKUnstructuredData."

    return is_valid_unstructured_mesh(dataset, repeat_cells)
end

function is_valid(dataset::VTKStructuredData)
    T = "VTKStructuredData"
    _dim = dim(dataset)
    in(_dim, [2,3]) && _dim == length(size(dataset.point_coords)[2:end]) || return false, "$T: Only 2 and 3 dimensional mesh data are allowed. First index iterates over the coordinate components, and the remaining indices iterate over the extents of the mesh."

    _extents = extents(dataset)
    all([_extents[i] > 0 for i in 1:length(_extents)]) || return false, "$T: Extent must be at least 1 in each dimensions."

    _keys = collect(keys(dataset.point_data))
    length(_keys) == 0 || begin
        _out = true
        for i in _keys
            var_extents = size(dataset.point_data[i])
            (_out = var_extents == _extents) || (_out = var_extents[2:end] == _extents) || break
        end
        _out
    end || return false, "$T: The point extents of the mesh and one of the point data fields are not consistent."

    cell_extents = (([_extents...] .- 1)...)
    _keys = collect(keys(dataset.cell_data))
    length(_keys) == 0 || begin
        _out = true
        for i in _keys
            var_extents = size(dataset.cell_data[i])
            (_out = var_extents == cell_extents) || (_out = var_extents[2:end] == cell_extents) || break
        end
        _out
    end || return false, "$T: The cell extents of the mesh and one of the cell data fields are not consistent."

    return true, ""
end

function is_valid{S}(dataset::VTKRectilinearData{S})
    T = "VTKRectilinearData"
    in(length(dataset.point_coords), [2,3]) || return false, "$T: Only 2 and 3 dimensional mesh data are allowed."

    !(S <: Real) || begin
        _out = true
        for i in length(dataset.point_coords), j in 1:length(dataset.point_coords[i])-1
            (_out = dataset.point_coords[i][j] < dataset.point_coords[i][j+1]) || break 
        end
        _out
    end || return false, "$T: Coordinate vectors should be in ascending order."

    _extents = extents(dataset)
    all([_extents[i] > 0 for i in 1:length(_extents)]) || return false, "$T: Extent must be at least 1 in each dimensions."
    _keys = collect(keys(dataset.point_data))
    length(_keys) == 0 || begin
        _out = true
        for i in _keys
            var_extents = size(dataset.point_data[i])
            (_out = var_extents == _extents) || (_out = var_extents[2:end] == _extents) || break
        end
        _out
    end || return false, "$T: The point extents of the mesh and one of the point data fields are not consistent."

    _cell_extents = cell_extents(dataset)
    _keys = collect(keys(dataset.cell_data))
    length(_keys) == 0 || begin
        _out = true
        for i in _keys
            var_extents = size(dataset.cell_data[i])
            (_out = var_extents == _cell_extents) || (_out = var_extents[2:end] == _cell_extents) || break
        end
        _out
    end || return false, "$T: The cell extents of the mesh and one of the cell data fields are not consistent."

    return true, ""
end

function is_valid{S}(dataset::VTKUniformRectilinearData{S})
    T = "VTKUniformRectilinearData"
    !(S <: Real) || begin 
        _out = true
        for i in 1:length(dataset.spacing)
            (_out = dataset.spacing[i] > 0) || return false, "Spacing must be positive for all dimensions."
        end
        _out
    end
    _extents = extents(dataset)
    all([_extents[i] > 0 for i in 1:length(_extents)]) || return false, "$T: Extent must be at least 1 in each dimensions."

    in(length(_extents), [2,3]) || return false, "$T: Only 2 and 3 dimensional mesh data are allowed."
    length(_extents) == length(dataset.spacing) == length(dataset.origin) || return false, "$T: The length of \"origin\", the length of \"spacing\", and the length of \"extents\" are not consistent."

    _keys = collect(keys(dataset.point_data))
    length(_keys) == 0 || begin
        _out = true
        for m in _keys
            var_extents = size(dataset.point_data[m])
            (_out = var_extents == _extents) || (_out = var_extents[2:end] == _extents) || break
        end
        _out
    end || return false, "$T: The point extents of the mesh and one of the point data fields are not consistent."

    cell_extents = (([_extents...] .- 1)...)
    _keys = collect(keys(dataset.cell_data))
    length(_keys) == 0 || begin
        _out = true
        for m in _keys
            var_extents = size(dataset.cell_data[m])
            (_out = var_extents == cell_extents) || (_out = var_extents[2:end] == cell_extents) || break
        end
        _out
    end || return false, "$T: The cell extents of the mesh and one of the cell data fields are not consistent."

    return true, ""
end

function is_valid(dataset::VTKMultiblockData) 
    for block in dataset
        valid, _error = is_valid(block)
        valid ? continue : return false, _error
    end
    return true, ""
end

function is_valid(dataset::VTKTimeSeriesData)
    T = "VTKTimeSeriesData"
    timemarkers = dataset.timemarkers
    data = dataset.data
    begin
        _out = true
        for i in 1:length(timemarkers)-1
            (_out = timemarkers[i] < timemarkers[i+1]) || break
        end
        _out
    end || return false, "$T: Time markers should be in ascending order."
    
    length(timemarkers) == length(data) || return false, "$T: The number of time markers and the length of data must be equal."

    for block in dataset
        valid, _error = is_valid(block)
        valid ? continue : return false, _error
    end
    return true, ""
end
