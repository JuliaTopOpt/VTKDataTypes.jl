__precompile__()
module VTKDataTypes

using WriteVTK
using PyGen
using GeometryTypes
using Colors
using Iterators

include("vtkcelltypes.jl")
include("types.jl")
include("utils.jl")
include("validation.jl")
include("typeconversion.jl")
include("decompose.jl")
include("triangulate.jl")
include("GLMesh.jl")
#include("surface.jl")

export  AbstractVTKData, 
        AbstractStaticVTKData, 
        AbstractTimeSeriesVTKData,
        AbstractVTKMultiblockData, 
        AbstractVTKSimpleData,
        AbstractVTKUnstructuredData, 
        AbstractVTKStructuredData, 
        AbstractVTKRectilinearData, 
        VTKUnstructuredData, 
        VTKPolyData, 
        VTKStructuredData, 
        VTKRectilinearData, 
        VTKUniformRectilinearData, 
        VTKImageData,
        VTKMultiblockData, 
        VTKTimeSeriesData, 
        promote_rule, 
        size, 
        length, 
        getindex, 
        endof, 
        start, 
        next, 
        done, 
        same_geometry, 
        ==, 
        same_ordered_geometry, 
        same_geometry_shape, 
        same_data_shape, 
        coherent, extents, 
        cell_extents, 
        dim, 
        num_of_points, 
        num_of_cells, 
        num_of_point_vars, 
        num_of_cell_vars, 
        cell_type, 
        cell_connectivity, 
        has_var, 
        var_dim, 
        is_homogeneous, 
        filter_cells!, 
        keep_volume_cells_only!, 
        keep_face_cells_only!, 
        get_cell_ids, 
        get_lowest_index, 
        get_highest_index, 
        is_valid_cell, 
        add_new_cell!, 
        remove_cell!, 
        add_point_id_offset!, 
        append!, 
        append, 
        num_of_blocks, 
        insert_new_block!, 
        remove_block!, 
        timespan, 
        num_of_timesteps, 
        insert_timed_data!, 
        remove_timed_data!, 
        is_valid, 
        convert, 
        simple_block_generator, 
        increase_resolution!,
        triangular, 
        triangulate, 
        decompose, 
        GLMesh, 
        bb, 
        pseudo_center
end
