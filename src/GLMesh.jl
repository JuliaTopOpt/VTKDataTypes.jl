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
    return 0.0, 0.0, 0.0, v
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

function GLMesh(dataset::VTKUnstructuredData; color::String="", component::Int=-1, opacity::Float64=1.0, color_scheme = BlackScheme())
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])
    tri_cell_data = color in keys(dataset.cell_data)
    tri_dataset = triangulate(dataset, tri_cell_data, GLTriangle)
    tri_dataset = VTKDataTypes.duplicate_vertices(tri_dataset)
    if haskey(dataset.cell_data, color)
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

function GLMesh(dataset::VTKPolyData; color::String="", component::Int=1, opacity::Float64=1.0)
    filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])

    faces = GLTriangle[]
    for i in 1:length(dataset.cell_connectivity)
        append!(faces, triangulate_cell_glmesh(dataset.cell_connectivity[i], dataset.cell_types[i]))
    end

    if color != ""
        color_variable = dataset.point_data[color]
        _var_dim = var_dim(dataset, color, "Point")
        if component == -1 || component == 1 && _var_dim == 1
            if _var_dim == 1
                _values = color_variable
                cmin = minimum(_values)
                cmax = maximum(_values)
            else
                @views _values = [norm(color_variable[:,i]) for i in 1:num_of_points(dataset)]
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
        if all(x -> isapprox(_values[1], x, atol = eps(Float32)), _values)
            scaled_value = zeros(length(values))
        else
            scaled_value = (_values .- cmin) ./ (cmax - cmin)
        end
        colors = RGBA{Float32}[RGBA{Float32}(color_scheme(scaled_value[i])...)
            for i=1:num_of_points(dataset)]
    else
        colors = RGBA{Float32}[RGBA{Float32}(1.0,1.0,1.0)
            for i=1:num_of_points(dataset)]
    end

    @views vertices = [Point{3, T}(dataset.point_coords[:,i]...) 
        for i in 1:num_of_points(dataset)]

    return GLNormalVertexcolorMesh(vertices=vertices, faces=faces, color=colors)
end

function GLMesh(_dataset::AbstractVTKStructuredData; color::String="", component::Int=1, opacity::Float64=1.0, color_scheme=BlackScheme())
    dataset = VTKStructuredData(_dataset)
    pextents = extents(dataset)
    _dim = dim(dataset)

    if _dim == 2
        faces = decompose_to_glmesh_2d(dataset)
    elseif _dim == 3
        faces = decompose_to_glmesh_3d(dataset)
    else
        throw("Unsupported dimension.")
    end

    if color != ""
        _color_variable = dataset.point_data[color]
        _var_dim = var_dim(dataset, color, "Point")
        if _var_dim == 1
            color_variable = reshape(_color_variable, (num_of_points(dataset),))
        else
            color_variable = reshape(_color_variable, (_var_dim, num_of_points(dataset)))
        end

        if component == -1 || component == 1 && _var_dim == 1
            if _var_dim == 1
                _values = color_variable
                cmin = minimum(_values)
                cmax = maximum(_values)
            else
                @views _values = [norm(color_variable[:,i]) for i in 1:num_of_points(dataset)]
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
        if all(x -> isapprox(_values[1], x, atol = eps(Float32)), _values)
            scaled_value = ones(length(values))
        else
            scaled_value = (_values .- cmin) ./ (cmax - cmin)
        end
        colors = RGBA{Float32}[RGBA{Float32}(color_scheme(scaled_value[i])...)
            for i=1:num_of_points(dataset)]
    else
        colors = RGBA{Float32}[RGBA{Float32}(0.1,0.1,0.1)
            for i=1:num_of_points(dataset)]
    end

    @views vertices = vec([Point{3, Float32}(dataset.point_coords[:,cind...]...) 
        for cind in Iterators.product([1:pextents[i] for i in 1:length(pextents)]...)])

    return GLNormalVertexcolorMesh(vertices=vertices, faces=faces, color=colors)
end

function decompose_to_glmesh_2d(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    faces = GLTriangle[]
    for cind in Iterators.product(1:cextents[1], 1:cextents[2])
        quad_cc = cell_connectivity(dataset, cind)
        push!(faces, GLTriangle(quad_cc[1], quad_cc[2], quad_cc[3]))
        push!(faces, GLTriangle(quad_cc[1], quad_cc[3], quad_cc[4]))
    end
    return faces
end

function decompose_to_glmesh_3d(dataset::VTKStructuredData)
    cextents = cell_extents(dataset)
    faces = GLTriangle[]
    for cind in Iterators.product(1:cextents[1], 1:cextents[2], 1:cextents[3])
        hexa_cc = cell_connectivity(dataset, cind)
        push!(faces, GLTriangle(hexa_cc[1], hexa_cc[2], hexa_cc[3]))
        push!(faces, GLTriangle(hexa_cc[1], hexa_cc[3], hexa_cc[4]))

        push!(faces, GLTriangle(hexa_cc[1], hexa_cc[2], hexa_cc[6]))
        push!(faces, GLTriangle(hexa_cc[1], hexa_cc[6], hexa_cc[5]))

        push!(faces, GLTriangle(hexa_cc[4], hexa_cc[1], hexa_cc[5]))
        push!(faces, GLTriangle(hexa_cc[4], hexa_cc[5], hexa_cc[8]))

        if cind[1] == cextents[1]
            push!(faces, GLTriangle(hexa_cc[2], hexa_cc[3], hexa_cc[7]))
            push!(faces, GLTriangle(hexa_cc[2], hexa_cc[7], hexa_cc[6]))
        end
        if cind[2] == cextents[2]
            push!(faces, GLTriangle(hexa_cc[4], hexa_cc[3], hexa_cc[7]))
            push!(faces, GLTriangle(hexa_cc[4], hexa_cc[7], hexa_cc[8]))
        end
        if cind[3] == cextents[3]
            push!(faces, GLTriangle(hexa_cc[5], hexa_cc[6], hexa_cc[7]))
            push!(faces, GLTriangle(hexa_cc[5], hexa_cc[7], hexa_cc[8]))
        end
    end

    return faces
end
