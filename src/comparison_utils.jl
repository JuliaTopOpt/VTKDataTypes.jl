function same_geometry(a::AbstractVTKUnstructuredData, b::AbstractVTKUnstructuredData)
    dim(a) != dim(b) && return false
    num_of_points(a) != num_of_points(b) && return false
    num_of_cells(a) != num_of_cells(b) && return false

    _num_of_points = num_of_points(a)
    b_image_point_inds = zeros(Int, _num_of_points)
    for i in 1:_num_of_points
        for j in 1:_num_of_points
            @views if b.point_coords[:, i] == a.point_coords[:, j]
                b_image_point_inds[i] = j
                break
            end
        end
        if b_image_point_inds[i] == 0
            return false
        end
    end

    function map_cell_connectivity(b_cell_connectivity)
        b_cell_connectivity_mapped = zeros(b_cell_connectivity)
        for i in 1:length(b_cell_connectivity)
            b_cell_connectivity_mapped[i] = b_image_point_inds[b_cell_connectivity[i]]
        end
        return b_cell_connectivity_mapped
    end

    _num_of_cells = num_of_cells(a)
    b_cell_connectivity_mapped = map(map_cell_connectivity, b.cell_connectivity)
    b_image_cell_inds = zeros(Int, _num_of_cells)
    for i in 1:_num_of_cells
        for j in 1:_num_of_cells
            if similar_cells(
                b_cell_connectivity_mapped[i],
                cell_type(b, i),
                a.cell_connectivity[j],
                cell_type(a, j),
            )
                b_image_cell_inds[i] = j
                break
            end
        end
        if b_image_cell_inds[i] == 0
            return false
        end
    end
    return true
end
function same_geometry(a::VTKStructuredData, b::VTKStructuredData)
    return a.point_coords == b.point_coords
end
function same_geometry(a::VTKRectilinearData, b::VTKRectilinearData)
    return extents(a) == extents(b) && all(a.point_coords .== b.point_coords)
end
function same_geometry(a::VTKUniformRectilinearData, b::VTKUniformRectilinearData)
    return a.spacing == b.spacing && extents(a) == extents(b)
end
function same_geometry(a::AbstractStaticVTKData, b::AbstractStaticVTKData)
    return same_geometry(promote(a, b)...)
end
function same_ordered_geometry(
    a::AbstractVTKUnstructuredData, b::AbstractVTKUnstructuredData
)
    dim(a) != dim(b) && return false
    num_of_points(a) != num_of_points(b) && return false
    num_of_cells(a) != num_of_cells(b) && return false

    _num_of_points = num_of_points(a)
    @views for i in 1:_num_of_points
        b.point_coords[:, i] == a.point_coords[:, i] || return false
    end

    _num_of_cells = num_of_cells(a)
    for i in 1:_num_of_cells
        similar_cells(
            a.cell_connectivity[i], cell_type(a, i), b.cell_connectivity[i], cell_type(b, i)
        ) || return false
    end
    return true
end
function same_ordered_geometry(a::AbstractVTKMultiblockData, b::AbstractVTKMultiblockData)
    length(a) == length(b) || return false
    for i in 1:length(a)
        typeof(a[i]) == typeof(b[i]) && (
            if typeof(a[i]) <: AbstractVTKMultiblockData ||
               typeof(a[i]) <: AbstractVTKUnstructuredData
                same_ordered_geometry(a[i], b[i])
            else
                a[i] == b[i]
            end
        ) || return false
    end
    return true
end

function same_geometry_shape(a::VTKStructuredData, b::VTKStructuredData)
    return size(a.point_coords) == size(b.point_coords)
end
function same_geometry_shape(a::VTKRectilinearData, b::VTKRectilinearData)
    return extents(a) == extents(b)
end
function same_geometry_shape(a::VTKUniformRectilinearData, b::VTKUniformRectilinearData)
    return extents(a) == extents(b)
end
function same_geometry_shape(a::AbstractVTKUnstructuredData, b::AbstractVTKUnstructuredData)
    dim(a) != dim(b) && return false
    num_of_points(a) != num_of_points(b) && return false
    num_of_cells(a) != num_of_cells(b) && return false

    _num_of_cells = num_of_cells(a)
    for i in 1:_num_of_cells
        similar_cells(
            a.cell_connectivity[i], cell_type(a, i), b.cell_connectivity[i], cell_type(b, i)
        ) || return false
    end
    return true
end
same_geometry_shape(a::AbstractStaticVTKData, b::AbstractStaticVTKData) = false

function same_ordered_geometry_shape(
    a::AbstractVTKMultiblockData, b::AbstractVTKMultiblockData
)
    length(a) == length(b) || return false
    for i in 1:length(a)
        typeof(a[i]) == typeof(b[i]) || return false
        if typeof(a[i]) <: AbstractVTKMultiblockData
            same_ordered_geometry_shape(a[i], b[i]) || return false
        else
            same_geometry_shape(a[i], b[i]) || return false
        end
    end
    return true
end

function same_data_shape(a::T, b::T) where {T<:AbstractVTKSimpleData}
    a_keys = keys(a.point_data)
    b_keys = keys(b.point_data)
    for i in a_keys
        in(i, b_keys) && size(a.point_data[i]) == size(b.point_data[i]) || return false
    end
    a_keys = keys(a.cell_data)
    b_keys = keys(b.cell_data)
    for i in a_keys
        in(i, b_keys) && size(a.cell_data[i]) == size(b.cell_data[i]) || return false
    end
    return true
end
same_data_shape(a::AbstractVTKStructuredData, b::AbstractVTKUnstructuredData) = false
same_data_shape(a::AbstractVTKMultiblockData, b::AbstractVTKSimpleData) = false

function same_data_shape(a::T, b::T) where {T<:AbstractVTKMultiblockData}
    length(a) == length(b)
    for i in 1:length(a)
        same_data_shape(a[i], b[i]) || return false
    end
    return true
end
