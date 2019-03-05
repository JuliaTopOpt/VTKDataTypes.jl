struct ColouredScheme{T}
    opacity::T
end
function (scheme::ColouredScheme)(v)
    0.0 <= v <= 1.0 || println("Warning: color value must be scaled between 0 and 1.")
    if 0.0 <= v < 0.25
        b = 1.0
        g = v/0.25
        r = 0.0
    elseif 0.25 <= v < 0.5
        b = 1.0 - (v - 0.25)/0.25
        g = 1.0
        r = 0.0
    elseif 0.5 <= v < 0.75
        b = 0.0
        g = 1.0
        r = (v - 0.5)/0.25
    else
        b = 0.0
        g = 1.0 - (v - 0.75)/0.25
        r = 1.0
    end
    return r, g, b, scheme.opacity
end

struct BlackScheme end
function (::BlackScheme)(v)
    0.0 <= v <= 1.0 || println("Warning: color value must be scaled between 0 and 1.")
    return 0.7, 0.7, 0.7, v
end

function get_jl_mapped_colors(dataset::T, color_variable_name, component, opacity, color_scheme = BlackScheme(), scale = true) where {T<:AbstractVTKSimpleData}
    if color_variable_name in keys(dataset.point_data)
        point_color = true
        _size = num_of_points(dataset)
        _color_variable = dataset.point_data[color_variable_name]
        _var_dim = var_dim(dataset, color_variable_name, "Point")
    elseif color_variable_name in keys(dataset.cell_data)
        point_color = false
        _size = num_of_cells(dataset)
        _color_variable = dataset.cell_data[color_variable_name]
        _var_dim = var_dim(dataset, color_variable_name, "Cell")
    else
        _size = num_of_points(dataset)
        _color_variable = nothing
    end

    if T <: AbstractVTKStructuredData
        if _var_dim == 1
            color_variable = reshape(_color_variable, (_size,))
        else
            color_variable = reshape(_color_variable, (_var_dim, _size))
        end
    else
        color_variable = _color_variable
    end
    
    if color_variable !== nothing
        if component == -1 || component == 1 && _var_dim == 1
            if _var_dim == 1
                _values = color_variable
                cmin = minimum(_values)
                cmax = maximum(_values)
            else
                @views _values = [norm(color_variable[:,i]) for i in 1:_size]
                cmin = minimum(_values)
                cmax = maximum(_values)
            end
        else
            if 1 <= component <= _var_dim
                @views _values = color_variable[component,:]
                cmin = minimum(_values)
                cmax = maximum(_values)
            else
                throw("Cannot use component $component of a $color. $color only has $_var_dim components.")
            end
        end
    else
        cmin = cmax = 1.0
        _values = ones(_size)
    end
    jl_mapped_colors = Vector{RGBA{Float32}}(undef, _size)
    if !(isapprox(cmax, cmin, atol = eps(Float32)))
        scaled_values = (_values .- cmin) ./ (cmax - cmin)
    else
        scaled_values = (similar(_values) .= clamp(cmin, 0, 1))
    end
    function fill_jl_mapped_colors(i)
        jl_mapped_colors[i] = RGBA{Float32}((color_scheme(scaled_values[i]))...)
    end
    map(fill_jl_mapped_colors, 1:_size)
    return jl_mapped_colors, cmin, cmax
end

function GLMesh(dataset::AbstractVTKUnstructuredData; color::String="", component::Int=-1, opacity::Float64=1.0, color_scheme = BlackScheme())
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])
    tri_cell_data = color in keys(dataset.cell_data)
    tri_dataset = triangulate(dataset, tri_cell_data, GLTriangle)
    tri_dataset = VTKDataTypes.duplicate_vertices(tri_dataset)
    if haskey(tri_dataset.cell_data, color)
        celldata_to_pointdata!(tri_dataset)
    end
    color_vec, cmin, cmax = get_jl_mapped_colors(   tri_dataset, 
                                                    color, 
                                                    component, 
                                                    opacity, 
                                                    color_scheme
                                                )
    T = eltype(dataset.point_coords)
    vertices = map(1:num_of_points(tri_dataset)) do i
        @views coord = tri_dataset.point_coords[:,i]
        if length(coord) == 2
            Point{3, T}(coord..., zero(T))
        else
            Point{3, T}(coord...)
        end
    end
    faces = tri_dataset.cell_connectivity
    return GLNormalVertexcolorMesh(vertices=vertices, faces=faces), color_vec
end
