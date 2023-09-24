abstract type AbstractVTKData end

abstract type AbstractStaticVTKData <: AbstractVTKData end
abstract type AbstractTimeSeriesVTKData{TTime, TData} <: AbstractVTKData end

abstract type AbstractVTKSimpleData <: AbstractStaticVTKData end
#abstract AbstractVTKAMRData <: AbstractStaticVTKData
abstract type AbstractVTKMultiblockData <: AbstractStaticVTKData end

abstract type AbstractVTKUnstructuredData <: AbstractVTKSimpleData end
abstract type AbstractVTKStructuredData <: AbstractVTKSimpleData end
abstract type AbstractVTKRectilinearData <: AbstractVTKStructuredData end

mutable struct VTKUnstructuredData{TCoord, TConnect, TPointData, TCellData} <: AbstractVTKUnstructuredData
    point_coords::Matrix{TCoord}
    cell_types::Vector{Int}
    cell_connectivity::Vector{TConnect}
    point_data::Dict{String, TPointData}
    cell_data::Dict{String, TCellData}
end
function VTKUnstructuredData(   point_coords, 
                                cell_types, 
                                cell_connectivity, 
                                point_data, 
                                cell_data, 
                                validate
                            )
    dataset = VTKUnstructuredData(  point_coords, 
                                    cell_types, 
                                    cell_connectivity, 
                                    point_data, 
                                    cell_data
                                )
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end

function VTKUnstructuredData(   point_coords, 
                                cell_types = Int[], 
                                cell_connectivity = Vector{Int}[], 
                                validate = false
                            )
    
    point_data = Dict{String, Array}()
    cell_data = Dict{String, Array}()    
    return VTKUnstructuredData( point_coords, 
                                cell_types, 
                                cell_connectivity, 
                                point_data, 
                                cell_data, 
                                validate
                              )
end

mutable struct VTKPolyData{TCoord, TConnect, TPointData, TCellData} <: AbstractVTKUnstructuredData
    point_coords::Matrix{TCoord}
    cell_types::Vector{Int}
    cell_connectivity::Vector{TConnect}
    point_data::Dict{String, TPointData}
    cell_data::Dict{String, TCellData}
end
function VTKPolyData(   point_coords, 
                        cell_types, 
                        cell_connectivity, 
                        point_data, 
                        cell_data, 
                        validate
                    )
    dataset = VTKPolyData(  point_coords, 
                            cell_types, 
                            cell_connectivity, 
                            point_data, 
                            cell_data
                        )
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end

function VTKPolyData(   point_coords, 
                        cell_types = Int[], 
                        cell_connectivity = Vector{Int}[], 
                        validate::Bool = false
                    )

    point_data = Dict{String, Array}()
    cell_data = Dict{String, Array}()
    return VTKPolyData( point_coords, 
                        cell_types, 
                        cell_connectivity, 
                        point_data,
                        cell_data,
                        validate
                      )
end

mutable struct VTKStructuredData{TCoord, TArray <: Array{TCoord}, TPointData, TCellData} <: AbstractVTKStructuredData
    point_coords::TArray
    point_data::Dict{String, TPointData}
    cell_data::Dict{String, TCellData}
end
function VTKStructuredData( point_coords, 
                            point_data, 
                            cell_data, 
                            validate::Bool
                          )
    dataset = VTKStructuredData(point_coords, point_data, cell_data)
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end
function VTKStructuredData(point_coords, validate::Bool)
    point_data = Dict{String, Array}()
    cell_data = Dict{String, Array}()
    return VTKStructuredData(point_coords, point_data, cell_data, validate)
end

mutable struct VTKRectilinearData{TCoord, N, TPointData, TCellData} <: AbstractVTKRectilinearData
    point_coords::NTuple{N, Vector{TCoord}}
    point_data::Dict{String, TPointData}
    cell_data::Dict{String, TCellData}
end
function VTKRectilinearData(point_coords, point_data, cell_data, validate)
    dataset = VTKRectilinearData(point_coords, point_data, cell_data)
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end
function VTKRectilinearData(point_coords, validate::Bool = false)
    point_data = Dict{String, Array}()
    cell_data = Dict{String, Array}()    
    VTKRectilinearData(point_coords, point_data, cell_data, validate)
end

mutable struct VTKUniformRectilinearData{TCoord, N, TPointData, TCellData} <: AbstractVTKRectilinearData
    origin::NTuple{N, TCoord}
    spacing::NTuple{N, TCoord}
    extents::NTuple{N, Int}
    point_data::Dict{String, TPointData}
    cell_data::Dict{String, TCellData}
end
function VTKUniformRectilinearData( origin, 
                                    spacing, 
                                    extents, 
                                    point_data, 
                                    cell_data, 
                                    validate::Bool
                                  )
    dataset = VTKUniformRectilinearData(origin, 
                                        spacing, 
                                        extents, 
                                        point_data, 
                                        cell_data
                                       )
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end
function VTKUniformRectilinearData( origin, 
                                    spacing, 
                                    extents, 
                                    validate::Bool = false,
                                  )
    point_data = Dict{String, Array}()
    cell_data = Dict{String, Array}()
    return VTKUniformRectilinearData(   origin, 
                                        spacing, 
                                        extents, 
                                        point_data, 
                                        cell_data, 
                                        validate
                                    )
end
const VTKImageData = VTKUniformRectilinearData

#=
type VTKBergerOligerAMRData <: AbstractVTKAMRData
    data::VTKUniformRectilinearData
    children::Vector{VTKBergerOligerAMRData}
    children_boundary_cells::Matrix{Int}
end
=#

mutable struct VTKMultiblockData{N, TBlocks <: Tuple{Vararg{AbstractStaticVTKData, N}}} <: AbstractVTKMultiblockData
    blocks::TBlocks
end

mutable struct VTKTimeSeriesData{TTime, TData <: AbstractStaticVTKData} <: AbstractTimeSeriesVTKData{TTime, TData}
    timemarkers::Vector{TTime}
    data::Vector{TData}
end
function VTKTimeSeriesData(timemarkers, data, validate)
    dataset = VTKTimeSeriesData(timemarkers, data)
    if validate
        valid, _error = is_valid(dataset)
        valid ? (return dataset) : throw(_error)
    end
    return dataset
end
