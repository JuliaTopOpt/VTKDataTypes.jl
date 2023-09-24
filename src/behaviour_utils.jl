function filter_cells!(
    dataset::AbstractVTKUnstructuredData, cell_types_to_remove::Vector{Int}
)
    inds = findall(i -> i in cell_types_to_remove, dataset.cell_types)
    deleteat!(dataset.cell_types, inds)
    deleteat!(dataset.cell_connectivity, inds)
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            deleteat!(dataset.cell_data[m], inds)
        else
            dataset.cell_data[m] = dataset.cell_data[m][
                :, setdiff(1:num_of_cells(dataset), inds)
            ]
        end
    end
    return dataset
end

function keep_volume_cells_only!(dataset::VTKUnstructuredData)
    return filter_cells!(dataset, [POINT_CELLS; LINE_CELLS; FACE_CELLS])
end
function keep_face_cells_only!(dataset::VTKPolyData)
    return filter_cells!(dataset, [POINT_CELLS; LINE_CELLS])
end
function add_new_cell!(
    dataset::AbstractVTKUnstructuredData,
    _cell_type::Int,
    _cell_connectivity,
    _filter::Bool=true,
)
    TData = typeof(dataset)
    begin
            out = true
            for i in 1:length(dataset.cell_connectivity)
                (
                    out =
                        !similar_cells(
                            dataset.cell_connectivity[i],
                            cell_type(dataset, i),
                            _cell_connectivity,
                            _cell_type,
                        )
                ) || break
            end
            out
        end ||
        _filter && (return num_of_cells(dataset)) ||
        throw("$(TData): Repeat cells are not allowed.")

    push!(dataset.cell_connectivity, _cell_connectivity)
    push!(dataset.cell_types, _cell_type)
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            dataset.cell_data[m] = [dataset.cell_data[m]; [0]]
        else
            dataset.cell_data[m] =
                [dataset.cell_data[m] zeros(size(dataset.cell_data[m], 1))]
        end
    end
    return num_of_cells(dataset)
end

function remove_cell!(dataset::AbstractVTKUnstructuredData, cell_ind::Int)
    deleteat!(dataset.cell_connectivity, cell_ind)
    deleteat!(dataset.cell_types, cell_ind)
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            dataset.cell_data[m] = vcat(
                dataset.cell_data[m][1:(cell_ind - 1)],
                dataset.cell_data[m][(cell_ind + 1):end],
            )
        else
            dataset.cell_data[m] = hcat(
                dataset.cell_data[m][:, 1:(cell_ind - 1)],
                dataset.cell_data[m][:, (cell_ind + 1):end],
            )
        end
    end
    return num_of_cells(dataset)
end

function add_point_id_offset!(cell_connectivity, offset::Int)
    for i in 1:length(cell_connectivity)
        cell_connectivity[i] = cell_connectivity[i] .+ offset
    end
    return cell_connectivity
end

function append(datasets::AbstractVTKUnstructuredData...)
    dataset1 = deepcopy(datasets[1])
    for i in 2:length(datasets)
        append!(dataset1, datasets[i])
    end
    return dataset1
end

insert_new_block!(dataset::AbstractVTKMultiblockData, block) = push!(dataset.blocks, block)
remove_block!(dataset::AbstractVTKMultiblockData, ind::Int) = deleteat!(dataset.blocks, ind)
function remove_block!(dataset::AbstractVTKMultiblockData, block::AbstractVTKMultiblockData)
    return deleteat!(dataset.blocks, findin(dataset.blocks, block))
end
function insert_timed_data!(
    dataset::VTKTimeSeriesData{TTime,TData}, time::TTime, data::TData
) where {TTime<:Real,TData<:AbstractStaticVTKData}
    loc = searchsorted(dataset.timemarkers, time)
    if loc.start == length(dataset) + 1 && loc.stop == length(dataset)
        push!(dataset.timemarkers, time)
        push!(dataset.data, data)
    elseif loc.start == 1 && loc.stop == 0
        unshift!(dataset.timemarkers, time)
        unshift!(dataset.data, data)
    else
        insert!(dataset.timemarkers, loc.start, time)
        insert!(dataset.data, loc.start, data)
    end
end

function remove_timed_data!(
    dataset::VTKTimeSeriesData{TTime,TData}, _time::TTime
) where {TTime<:Real,TData<:AbstractStaticVTKData}
    ind = findin(dataset.timemarkers, _time)
    if typeof(ind) == Int
        deleteat!(dataset.timemarkers, ind)
        deleteat!(dataset.data, ind)
    else
        throw("No data at time $_time.")
    end
end

function remove_timed_data!(dataset::VTKTimeSeriesData, ind::Int)
    if ind >= 1 && ind <= length(dataset)
        deleteat!(dataset.timemarkers, ind)
        deleteat!(dataset.data, ind)
    else
        throw("Index out of bound!")
    end
end

function increase_resolution!(dataset::VTKTimeSeriesData, resolution::Int)
    _times = deepcopy(
        collect(zip(dataset.timemarkers[1:(end - 1)], dataset.timemarkers[2:end]))
    )
    for (t1, t2) in _times
        _step = (t2 - t1) / resolution
        for j in 1:(resolution - 1)
            VTKDataTypes.insert_timed_data!(
                dataset, Float64(t1 + j * _step), dataset[t1 + j * _step]
            )
        end
    end
    return nothing
end

@resumable function simple_block_generator(multiblock::VTKMultiblockData)
    for block in multiblock
        if typeof(block) <: VTKMultiblockData
            for b in simple_block_generator(block)
                @yield b
            end
        else
            @yield block
        end
    end
end

function dim3!(a::AbstractVTKUnstructuredData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        a.point_coords = [a.point_coords; zeros(1, num_of_points(a))]
        return nothing
    else
        throw("Invalid dimension.")
    end
end

function dim3!(a::VTKUniformRectilinearData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        a.origin = [a.origin; 0]
        a.spacing = [a.spacing; 0]
        a.extents = [a.extents; 1]
        for m in keys(a.point_data)
            a.point_data[m] = reshape(a.point_data[m], (size(a.point_data[m])..., 1))
        end
        for m in keys(a.cell_data)
            a.cell_data[m] = reshape(a.cell_data[m], (size(a.cell_data[m])..., 1))
        end

        return nothing
    else
        throw("Invalid dimension.")
    end
end

function dim3!(a::VTKRectilinearData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        a.point_coords = [a.point_coords; [0]]
        for m in keys(a.point_data)
            a.point_data[m] = reshape(a.point_data[m], (size(a.point_data[m])..., 1))
        end
        for m in keys(a.cell_data)
            a.cell_data[m] = reshape(a.cell_data[m], (size(a.cell_data[m])..., 1))
        end

        return nothing
    else
        throw("Invalid dimension.")
    end
end

function dim3!(a::VTKStructuredData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        a.point_coords = reshape(
            cat(1, a.point_coords, zeros(1, extents(a)...)), (3, extents(a), 1)
        )
        for m in keys(a.point_data)
            a.point_data[m] = reshape(a.point_data[m], (size(a.point_data[m])..., 1))
        end
        for m in keys(a.cell_data)
            a.cell_data[m] = reshape(a.cell_data[m], (size(a.cell_data[m])..., 1))
        end

        return nothing
    else
        throw("Invalid dimension.")
    end
end

function dim3!(a::VTKMultiblockData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        for block in simple_block_generator(a)
            dim3!(block)
        end
        return nothing
    else
        throw("Invalid dimension.")
    end
end

function dim3!(a::VTKTimeSeriesData)
    if dim(a) == 3
        return nothing
    elseif dim(a) == 2
        for timedblock in a
            dim3!(timedblock)
        end
        return nothing
    else
        throw("Invalid dimension.")
    end
end

function celldata_to_pointdata!(dataset::AbstractVTKUnstructuredData)
    point_coords = dataset.point_coords
    point_data = dataset.point_data

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            point_data[m] = zeros(num_of_points(dataset))
        else
            point_data[m] = zeros(_var_dim, num_of_points(dataset))
        end
    end

    point_counter = zeros(num_of_points(dataset))
    for i in 1:num_of_cells(dataset), j in dataset.cell_connectivity[i]
        point_counter[j] += 1
        for m in keys(dataset.cell_data)
            _var_dim = var_dim(dataset, m, "Cell")
            if _var_dim == 1
                point_data[m][j] += dataset.cell_data[m][i]
            else
                @views point_data[m][:, j] += dataset.cell_data[m][:, i]
            end
        end
    end

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            point_data[m] = point_data[m] ./ max.(point_counter, 1)
        else
            point_data[m] = point_data[m] ./ max.(point_counter', 1)
        end
    end
    empty!(dataset.cell_data)

    return nothing
end
