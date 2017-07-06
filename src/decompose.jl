include("decompose_cell.jl")
include("decompose_unstruct.jl")
include("decompose_struct.jl")

decompose(dataset::AbstractVTKMultiblockData, target::String="Faces") = decompose(VTKUnstructuredData(dataset), target)

function decompose(dataset::AbstractTimeSeriesVTKData, target::String="Faces")
    return VTKTimeSeriesData(dataset.timemarkers, [decompose(dataset[i], target) for 
        i in 1:length(dataset)])
end

