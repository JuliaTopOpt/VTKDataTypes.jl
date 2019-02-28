# VTKDataTypes
## Overview

VTKDataTypes.jl presents a Julia type system for representing and manipulating VTK data natively in Julia. **VTKDataTypes.jl only supports Julia v1.**

## Summary of capabilities

### VTK data representation

You can use VTKDataTypes.jl to create your own VTK data object of any of the following types: `VTKUnstructuredData`, `VTKPolyData`, `VTKStructuredData`, `VTKRectilinearData`, `VTKImageData`, `VTKMultiblockData`, and `VTKTimeSeriesData`. 2D and 3D data are both supported and conversion functions are defined where they make sense.

### Decomposing and triangulation

You can use the `decompose` function to decompose your VTK data object to faces or lines averaging cell data of all participating cells in the face or edge. The `triangulate` function is also defined for interfacing with other rendering and visualization modules which require a triangular mesh.

### Interfacing with GeometryTypes.jl

You can use the `GLMesh` function to change a VTK data object to a triangulated `GLNormalVertexcolorMesh` using one of the point data variable names to define the vertex color specified by the `color` option. Other options include `opacity`, and `component` to specify a vector data component, as opposed to the magnitude which is used by default. 

## VTK Data Types

### VTK Cell Types

The following is a record of the ID associated with the most common cell types as per VTK's convention:

1	: 	VTK_VERTEX

2	: 	VTK_POLY_VERTEX

3	:	VTK_LINE

4	:	VTK_POLY_LINE

5	:	VTK_TRIANGLE

6	:	VTK_TRIANGLE_STRIP

7	:	VTK_POLYGON

8	:	VTK_PIXEL

9	:	VTK_QUAD

10	:	VTK_TETRA 

11	:	VTK_VOXEL

12	:	VTK_HEXAHEDRON

13	:	VTK_WEDGE 

14	:	VTK_PYRAMID 

15	:	VTK_PENTAGONAL_PRISM

16	:	VTK_HEXAGONAL_PRISM 

For more, you can refer to `src/vtkcelltypes.jl`.

### VTKUnstructuredData{T}
`point_coords` : this is a `Matrix{T}` of size `(dims, points)` where `dims` is the number of dimensions and `points` is the number of points.

`cell_types` : this is a `Vector{Int}` holding the cell types using VTK's convention for cell ids.

`cell_connectivity` : this is a `Vector{Vector{Int}}` holding the cell connectivity of each cell as per VTK's convention for each cell type.

`point_data` : this is a `Dict{String, Array}` holding point-specific scalar and vector data. Scalar data arrays are `Vector{T}` and vector data arrays are `Matrix{T}` of size `(var_dims, points)` where `var_dims` is the number of dimensions of the vector variable.

`cell_data` : this is a `Dict{String, Array}` holding cell-specific scalar and vector data. Scalar data arrays are `Vector{T}` and vector data arrays are `Matrix{T}` of size `(var_dims, cells)` where `var_dims` is the number of dimensions of the vector variable and `cells` is the number of cells.

### VTKPolyData{T}

Has the same fields as `VTKUnstructuredData` but all cells are assumed to be of types: 1, 2, 3, 4, 5, 6, 7, 8, and/or 9.

### VTKStructuredData{T}

`point_coords` : this is an `Array{T,dims+1}`. The size of the array is `(dims, extents(dataset)...)` where `extents` is a function that returns the number of grid markers along each dimension.

`point_data` : this is a `Dict{String, Array}` holding point-specific scalar and vector data. Scalar data arrays are stored in `Array{T,dims}`. The size of a scalar data array is `extents(dataset)`. Vector data arrays are stored in `Array{T,dims+1}` of size `(var_dims, extents(dataset)...)`.

`cell_data` : this is a `Dict{String, Array}` holding cell-specific scalar and vector data. Scalar data arrays are stored in `Array{T,dims}`. The size of a scalar data array is `cell_extents(dataset)`, where `cell_extents` is a function that returns `extents(dataset) .- 1` resembling the number of cells along each dimension. Vector data arrays are stored in `Array{T,dims+1}` of size `(var_dims, cell_extents(dataset)...)`.

All cells are of type 9 (Quad) in 2D and 12 (Hexa) in 3D.

### VTKRectilinearData{T}

`point_coords` : this is a `NTuple{dim, Vector{T}}` holding the x, y (and z) coordinates of the `dim` dimensional rectilinear grid.

`point_data` : similar to `VTKStructuredData{T}`

`cell_data` : similar to `VTKStructuredData{T}`

All cells are of type 9 (Quad) in 2D and 12 (Hexa) in 3D.

### VTKUniformRectilinearData{T} aka VTKImageData{T}

`origin` : an `NTuple{dim, T}` that refers to the origin of the uniform rectilinear structure of dimension `dim`

`spacing` : an `NTuple{dim, T}` that refers to the spacing between grid lines along each axis

`extents` : an `NTuple{dim, Int}` that refers to the number of points along each axis

`point_data` : similar to `VTKStructuredData{T}`

`cell_data` : similar to `VTKStructuredData{T}`

All cells are of type 8 (Pixel) in 2D and 11 (Voxel) in 3D.

### VTKMultiblockData

`blocks` : a `Tuple{Vararg{AbstractStaticVTKData}}` which can store a mixture of all previous types as well as `VTKMultiblockData` recursively.

Iteration and indexing are defined for this type.

### VTKTimeSeriesData{TTime, TData <: AbstractStaticVTKData}

`timemarkers` : a `Vector{TTime}` of frame times.

`data` : a `Vector{TData}` of any one of the previous data types. The type must be consistent for all the time steps.

Iteration and indexing are defined for this type. Integer indexing will access the frames by their index in `data`. Frames can also be accessed by their time using floating point indexing. Automatic linear interpolation and constant extrapolation is done when indexing with time.

## Examples

Please refer to the `tests` folder for examples.
