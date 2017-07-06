abstract type AbstractVTKData end

abstract type AbstractStaticVTKData{T} <: AbstractVTKData end
abstract type AbstractTimeSeriesVTKData{S, T} <: AbstractVTKData end

abstract type AbstractVTKSimpleData{T} <: AbstractStaticVTKData{T} end
#abstract AbstractVTKAMRData <: AbstractStaticVTKData
abstract type AbstractVTKMultiblockData{T} <: AbstractStaticVTKData{T} end

abstract type AbstractVTKUnstructuredData{T} <: AbstractVTKSimpleData{T} end
abstract type AbstractVTKStructuredData{T} <: AbstractVTKSimpleData{T} end
abstract type AbstractVTKRectilinearData{T} <: AbstractVTKStructuredData{T} end

type VTKUnstructuredData{T} <: AbstractVTKUnstructuredData{T}
    point_coords::Matrix{T}
    cell_types::Vector{Int}
    cell_connectivity::Vector{Vector{Int}}
    point_data::Dict{String, Array{T}}
    cell_data::Dict{String, Array{T}}

    function VTKUnstructuredData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, validate=false) where T
        dataset = new(point_coords, cell_types, cell_connectivity, point_data, cell_data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end

function VTKUnstructuredData{T1, T2, T3}(point_coords::Matrix{T1}, cell_types::Vector{Int},
    cell_connectivity::Vector{Vector{Int}}, point_data::Dict{String, Array{T2}}, 
    cell_data::Dict{String, Array{T3}}, validate=false)
    
    T = promote_type(T1, T2, T3)
    if T <: Integer
        T = Float64
    end
    println("Here")
    VTKUnstructuredData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, validate)
end

function VTKUnstructuredData{_T}(point_coords::Matrix{_T}, cell_types::Vector{Int}=Int[],
    cell_connectivity::Vector{Vector{Int}}=[Int[]], validate=false)
    
    T = _T
    if T <: Integer
        T = Float64
    end
    point_data = Dict{String, Array{T}}()
    cell_data = Dict{String, Array{T}}()
    
    VTKUnstructuredData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, validate)
end

type VTKPolyData{T} <: AbstractVTKUnstructuredData{T}
    point_coords::Matrix{T}
    cell_types::Vector{Int}
    cell_connectivity::Vector{Vector{Int}}
    point_data::Dict{String, Array{T}}
    cell_data::Dict{String, Array{T}}

    function VTKPolyData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, validate=false) where T
        dataset = new(point_coords, cell_types, cell_connectivity, point_data, cell_data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end
function VTKPolyData{T1, T2, T3}(point_coords::Matrix{T1}, cell_types::Vector{Int},
    cell_connectivity::Vector{Vector{Int}}, point_data::Dict{String, Array{T2}}, cell_data::Dict{String, Array{T3}}, validate=false)
    
    T = promote_type(T1, T2, T3)
    if T <: Integer
        T = Float64
    end

    VTKPolyData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, validate)
end
function VTKPolyData{_T}(point_coords::Matrix{_T}, cell_types::Vector{Int}=Int[],
    cell_connectivity::Vector{Vector{Int}}=[Int[]], validate=false)

    T = _T
    if T <: Integer
        T = Float64
    end
    point_data = Dict{String, Array{T}}()
    cell_data = Dict{String, Array{T}}()
    VTKPolyData{T}(point_coords, cell_types, cell_connectivity, point_data, cell_data, dataset)
end

type VTKStructuredData{T} <: AbstractVTKStructuredData{T}
    point_coords::Array{T}
    point_data::Dict{String, Array{T}}
    cell_data::Dict{String, Array{T}}

    function VTKStructuredData{T}(point_coords, point_data, cell_data, validate=false) where T
        dataset = new(point_coords, point_data, cell_data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end
function VTKStructuredData{T1, N, T2, T3}(point_coords::Array{T1,N}, 
    point_data::Dict{String, Array{T2}}, cell_data::Dict{String, Array{T3}}, 
    validate=false)

    T = promote_type(T1, T2, T3)
    if T <: Integer
        T = Float64
    end

    VTKStructuredData{T}(point_coords, point_data, cell_data, validate)
end
function VTKStructuredData{_T}(point_coords::Array{_T}, validate=false)
    T = _T
    if T <: Integer
        T = Float64
    end
    point_data = Dict{String, Array{T}}()
    cell_data = Dict{String, Array{T}}()
    VTKStructuredData{T}(point_coords, point_data, cell_data, validate)
end

type VTKRectilinearData{T} <: AbstractVTKRectilinearData{T}
    point_coords::Vector{Vector{T}}
    point_data::Dict{String, Array{T}}
    cell_data::Dict{String, Array{T}}

    function VTKRectilinearData{T}(point_coords, point_data, cell_data, validate=false) where T
        dataset = new(point_coords, point_data, cell_data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end
function VTKRectilinearData{T1, T2, T3}(point_coords::Vector{Vector{T1}}, 
    point_data::Dict{String, Array{T2}}, cell_data::Dict{String, Array{T3}}, 
    validate=false)
    
    T = promote_type(T1, T2, T3)
    if T <: Integer
        T = Float64
    end

    VTKRectilinearData{T}(point_coords, point_data, cell_data, validate)
end
function VTKRectilinearData{_T}(point_coords::Vector{Vector{_T}}, validate=false)
    T = _T
    if T <: Integer
        T = Float64
    end
    point_data = Dict{String, Array{T}}()
    cell_data = Dict{String, Array{T}}()
    
    VTKRectilinearData{T}(point_coords, point_data, cell_data, validate)
end

type VTKUniformRectilinearData{T} <: AbstractVTKRectilinearData{T}
    origin::Vector{T}
    spacing::Vector{T}
    extents::Vector{Int}
    point_data::Dict{String, Array{T}}
    cell_data::Dict{String, Array{T}}

    function VTKUniformRectilinearData{T}(origin, spacing, extents, point_data, cell_data, validate=false) where T
        dataset = new(origin, spacing, extents, point_data, cell_data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end
function VTKUniformRectilinearData{T1, T2, T3, T4}(origin::Vector{T1}, 
    spacing::Vector{T2}, extents::Vector{Int}, point_data::Dict{String, Array{T3}}, 
    cell_data::Dict{String, Array{T4}}, validate=false)
    
    T = promote_type(T1, T2, T3, T4)
    if T <: Integer
        T = Float64
    end

    VTKUniformRectilinearData{T}(origin, spacing, extents, point_data, cell_data, validate)
end
function VTKUniformRectilinearData{T1, T2}(origin::Vector{T1}, spacing::Vector{T2}, extents::Vector{Int}, validate=false)
    T = promote_type(T1, T2)
    if T <: Integer
        T = Float64
    end
    point_data = Dict{String, Array{T}}()
    cell_data = Dict{String, Array{T}}()
    VTKUniformRectilinearData{T}(origin, spacing, extents, point_data, cell_data, validate)
end
const VTKImageData = VTKUniformRectilinearData

#=
type VTKBergerOligerAMRData <: AbstractVTKAMRData
    data::VTKUniformRectilinearData
    children::Vector{VTKBergerOligerAMRData}
    children_boundary_cells::Matrix{Int}
end
=#

type VTKMultiblockData{T} <: AbstractVTKMultiblockData{T}
    blocks::Vector{AbstractStaticVTKData{T}}
end

type VTKTimeSeriesData{S, T<:AbstractStaticVTKData} <: AbstractTimeSeriesVTKData{S, T}
    timemarkers::Vector{T}
    data::Vector{S}

    function VTKTimeSeriesData{S, T}(timemarkers, data, validate=false) where {S, T<:AbstractStaticVTKData}
        dataset = new(timemarkers, data)
        if validate
            valid, _error = is_valid(dataset)
            valid ? (return dataset) : throw(_error)
        end
        dataset
    end
end
VTKTimeSeriesData{T, S<:AbstractStaticVTKData}(timemarkers::Vector{T}, data::Vector{S}, validate=false) = 
    T <: Integer ? VTKTimeSeriesData{Float64, S}(timemarkers, data) : 
    VTKTimeSeriesData{T, S}(timemarkers, data, validate)
