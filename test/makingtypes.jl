# This is an example that shows you how to create VTK data types natively in Julia

using VTKDataTypes
using GeometryTypes: GeometryTypes

function create_image()
    origin = (0.0, 0.0, 0.0)
    spacing = (1.0, 2.0, 1.0)
    _extents = (2, 3, 5) #Must be Vector{Int}

    image = VTKImageData(origin, spacing, _extents)
    @test dim(image) == 3
    @test extents(image) == image.extents
    @test cell_extents(image) == image.extents .- 1

    image.point_data["Point scalar data"] = zeros(extents(image)...)
    image.cell_data["Cell scalar data"] = zeros(cell_extents(image)...)
    image.point_data["Point vector data"] = zeros(3, extents(image)...)
    image.cell_data["Cell vector data"] = zeros(3, cell_extents(image)...)

    @test cell_connectivity(image, (1, 1, 1)) == [1, 2, 3, 4, 7, 8, 9, 10]
    #Used to obtain the points which make up the cell in VTK convention
    @test cell_type(image, (1, 1, 1)) == 11
    #Used to obtain the VTK cell type

    @test image == deepcopy(image)
    return image
end

function create_rectilinear_data1()
    x = y = z = collect(range(-2; stop=2, length=5))

    rectilinear = VTKRectilinearData((x, y, z))
    @test dim(rectilinear) == 3
    @test extents(rectilinear) == (5, 5, 5)
    @test cell_extents(rectilinear) == (4, 4, 4)

    rectilinear.point_data["Point scalar data"] = zeros(extents(rectilinear)...)
    rectilinear.cell_data["Cell scalar data"] = zeros(cell_extents(rectilinear)...)
    rectilinear.point_data["Point vector data"] = zeros(3, extents(rectilinear)...)
    rectilinear.cell_data["Cell vector data"] = zeros(3, cell_extents(rectilinear)...)

    @test cell_connectivity(rectilinear, (1, 1, 1)) == [1, 2, 7, 6, 26, 27, 32, 31]
    #Used to obtain the points which make up the cell in VTK convention
    @test cell_type(rectilinear, (1, 1, 1)) == 12
    #Used to obtain the VTK cell type

    @test rectilinear == deepcopy(rectilinear)
    return rectilinear
end

function create_rectilinear_data2()
    image = create_image()
    rectilinear = VTKRectilinearData(image)

    @test dim(rectilinear) == dim(image) == 3
    @test extents(rectilinear) == extents(image) == (2, 3, 5)
    @test cell_extents(rectilinear) == cell_extents(image) == (1, 2, 4)

    return rectilinear
end

function create_structured_data()
    rectilinear = create_rectilinear_data1()
    structured = VTKStructuredData(rectilinear)
    @test size(structured.point_coords) == (3, extents(rectilinear)...) == (3, 5, 5, 5)

    @test extents(structured) == extents(rectilinear) == (5, 5, 5)
    @test cell_extents(structured) == cell_extents(rectilinear) == (4, 4, 4)
    @test dim(structured) == dim(rectilinear) == 3

    structured.point_data["Point scalar data"] = zeros(extents(structured)...)
    structured.cell_data["Cell scalar data"] = zeros(cell_extents(structured)...)
    structured.point_data["Point vector data"] = zeros(3, extents(structured)...)
    structured.cell_data["Cell vector data"] = zeros(3, cell_extents(structured)...)

    @test cell_connectivity(structured, (1, 1, 1)) == [1, 2, 7, 6, 26, 27, 32, 31]
    #Used to obtain the points which make up the cell in VTK convention
    @test cell_type(structured, (1, 1, 1)) == 12
    #Used to obtain the VTK cell type

    @test structured == deepcopy(structured)
    return structured
end

function create_unstructured_data()
    rectilinear = create_rectilinear_data1()
    structured = create_structured_data()
    unstruct1 = VTKUnstructuredData(rectilinear)
    unstruct2 = VTKUnstructuredData(structured)

    @test length(unstruct1.cell_connectivity) == length(unstruct1.cell_types)
    @test size(unstruct1.point_coords) == (dim(unstruct1), num_of_points(unstruct1))
    @test size(unstruct1.point_data["Point scalar data"], 1) ==
        size(unstruct1.point_coords, 2)
    @test size(unstruct1.cell_data["Cell vector data"]) == (3, num_of_cells(unstruct1))

    @test unstruct1 == unstruct2

    @test cell_connectivity(unstruct1, 1) == unstruct1.cell_connectivity[1]
    @test cell_type(unstruct1, 1) == unstruct1.cell_types[1]
    @test is_homogeneous(unstruct1)
    return unstruct1
end

function create_poly_data()
    rectilinear = create_rectilinear_data1()
    polydata = VTKPolyData(rectilinear)
    # automatically decomposes linear volume cells to faces and averages cell data

    @test length(polydata.cell_connectivity) == length(polydata.cell_types)
    @test size(polydata.point_coords) == (dim(polydata), num_of_points(polydata))
    @test size(polydata.point_data["Point scalar data"], 1) ==
        size(polydata.point_coords, 2)
    @test size(polydata.cell_data["Cell vector data"]) == (3, num_of_cells(polydata))

    @test polydata == polydata
    return polydata
end

function create_multiblock_data()
    mb = VTKMultiblockData((
        create_image(),
        create_rectilinear_data1(),
        create_structured_data(),
        create_unstructured_data(),
        create_poly_data(),
    ))
    @test length(mb) == 5
    for b in mb
        @test isa(b, AbstractVTKSimpleData)
    end
    mb2 = VTKMultiblockData((mb, create_image()))
    for b in simple_block_generator(mb2) #Recursively unwraps any VTKMultiblockData
        @test isa(b, AbstractVTKSimpleData)
    end

    return mb2[1]
end

function create_timeseries_data()
    timemarkers = [0, 0.5, 1, 1.5]

    #2D polydata
    x = y = [0.0, 1.0, 2.0]
    point_coords = (x, y)
    point_data = Dict{String,Array{Float64}}()
    cell_data = Dict{String,Array{Float64}}()
    a = VTKRectilinearData(point_coords, point_data, cell_data)
    a.cell_data["Color"] = reshape(rand(num_of_cells(a)), cell_extents(a))
    b = VTKPolyData(a)

    x = y = [3.0, 4.0, 5.0]
    point_coords = (x, y)
    c = VTKPolyData(VTKRectilinearData(point_coords, point_data, cell_data))
    c.cell_data["Color"] = rand(num_of_cells(c))

    m = VTKMultiblockData((b, c))
    timeseries = VTKTimeSeriesData([0.0], [m])

    no_of_timesteps = 5
    timestep = 0.5
    speed = 1
    i = 0
    while i < no_of_timesteps
        i += 1
        new_data = deepcopy(timeseries.data[i])

        #Random walk of first dataset
        new_data[1].point_coords =
            new_data[1].point_coords + vcat(
                fill(speed * rand() - speed / 2, (1, num_of_points(new_data[1]))),
                fill(speed * rand() - speed / 2, (1, num_of_points(new_data[1]))),
            )

        #Random walk of second dataset
        new_data[2].point_coords =
            new_data[2].point_coords + vcat(
                fill(speed * rand() - speed / 2, (1, num_of_points(new_data[2]))),
                fill(speed * rand() - speed / 2, (1, num_of_points(new_data[2]))),
            )

        insert_timed_data!(timeseries, i * timestep, new_data)
    end

    for ts in timeseries
        @test isa(ts, VTKMultiblockData)
    end

    @test isa(timeseries[1.1], VTKMultiblockData) #Indexing by time, does linear interpolation
    @test isa(timeseries[1], VTKMultiblockData) #Indexing by frame index

    return timeseries
end

function triangulation_and_decompose()
    rectilinear = create_rectilinear_data1()
    tridata = triangulate(rectilinear)
    @test triangular(tridata)
    @test isa(GLMesh(rectilinear), GeometryTypes.HomogenousMesh)
    @test isa(decompose(rectilinear, "Faces"), VTKPolyData)
    @test triangulate(decompose(rectilinear, "Faces")) == tridata

    return tridata
end

create_image()
create_rectilinear_data1()
create_rectilinear_data2()
create_structured_data()
create_unstructured_data()
create_poly_data()
create_multiblock_data()
create_timeseries_data()
triangulation_and_decompose()
