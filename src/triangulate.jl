include("triangulate_cell.jl")
include("triangulate_unstruct.jl")
include("triangulate_struct.jl")

triangulate(dataset::AbstractVTKMultiblockData) = triangulate(VTKUnstructuredData(dataset))

function triangulate(dataset::AbstractTimeSeriesVTKData)
    return VTKTimeSeriesData(
        dataset.timemarkers, [triangulate(dataset[i]) for i in 1:length(dataset)]
    )
end
