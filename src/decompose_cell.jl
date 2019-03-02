function decompose_cell(cell_connectivity, cell_type::Int; target::String="")
    if target == ""
        if in(cell_type, POINT_CELLS) || in(cell_type, LINE_CELLS)
            throw("Cell is already represented as points.")
        elseif in(cell_type, FACE_CELLS)
            target = "Lines"
        elseif in(cell_type, VOLUME_CELLS)
            target = "Faces"
        end
    end

    args = (cell_connectivity, target)
    if cell_type == 5
        return decompose_triangle(args...)
    elseif cell_type == 6
        return decompose_triangle_strip(args...)
    elseif cell_type == 8
        return decompose_pixel(args...)
    elseif cell_type == 9
        return decompose_quad(args...)
    elseif cell_type == 10
        return decompose_tetra(args...)
    elseif cell_type == 11
        return decompose_voxel(args...)
    elseif cell_type == 12
        return decompose_hexa(args...)
    elseif cell_type == 13
        return decompose_wedge(args...)
    elseif cell_type == 14
        return decompose_pyramid(args...)
    elseif cell_type == 22
        return decompose_quadratic_triangle(args...)
    elseif cell_type == 23
        return decompose_quadratic_quad(args...)
    elseif cell_type == 24
        return decompose_quadratic_tetra(args...)
    elseif cell_type == 25
        return decompose_quadratic_hexa(args...)
    else
        throw("Cell type $cell_type is not supported.")
    end
end

function decompose_triangle(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[cell_connectivity], Int[5]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ], 
            [ cell_connectivity[2], cell_connectivity[3] ], 
            [ cell_connectivity[3], cell_connectivity[1] ] ],
            fill(3,3)
    end
end

function decompose_triangle_strip(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[i], cell_connectivity[i+1], cell_connectivity[i+2] ] 
            for i in 1:length(cell_connectivity)-2], fill(5, length(cell_connectivity)-2)
    else
        return append!(
            Vector{Int}[ [ cell_connectivity[i], cell_connectivity[i+1] ]
                for i in 1:length(cell_connectivity)-1 ], 
            Vector{Int}[ [ cell_connectivity[i], cell_connectivity[i+2] ]
                for i in 1:length(cell_connectivity)-2]), 
            fill(3, 2*length(cell_connectivity)-3)
    end
end

function decompose_pixel(cell_connectivity, target::String)
    if target == "Faces"
        return [cell_connectivity], [8]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[4] ],
            [ cell_connectivity[4], cell_connectivity[3] ], 
            [ cell_connectivity[3], cell_connectivity[1] ] ], 
            fill(3, 4)
    end
end

function decompose_quad(cell_connectivity, target::String)
    if target == "Faces"
        return [cell_connectivity], [9]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ], 
            [ cell_connectivity[2], cell_connectivity[3] ],
            [ cell_connectivity[3], cell_connectivity[4] ],
            [ cell_connectivity[4], cell_connectivity[1] ] ],
            fill(3, 4)
    end
end

function decompose_tetra(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3] ], 
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[4] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[4] ],
            [ cell_connectivity[2], cell_connectivity[1], cell_connectivity[4] ] ], 
            fill(5, 4)
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[3] ],
            [ cell_connectivity[3], cell_connectivity[1] ],
            [ cell_connectivity[1], cell_connectivity[4] ], 
            [ cell_connectivity[2], cell_connectivity[4] ], 
            [ cell_connectivity[3], cell_connectivity[4] ] ], 
            fill(3, 6)
    end
end

function decompose_voxel(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[4], cell_connectivity[3] ], 
            [ cell_connectivity[5], cell_connectivity[6], cell_connectivity[8], cell_connectivity[7] ],
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[6], cell_connectivity[5] ], 
            [ cell_connectivity[2], cell_connectivity[4], cell_connectivity[8], cell_connectivity[6] ],
            [ cell_connectivity[4], cell_connectivity[3], cell_connectivity[7], cell_connectivity[8] ],
            [ cell_connectivity[3], cell_connectivity[1], cell_connectivity[5], cell_connectivity[7] ] ],
            fill(9, 6)
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[4] ],
            [ cell_connectivity[4], cell_connectivity[3] ], 
            [ cell_connectivity[3], cell_connectivity[1] ],
            [ cell_connectivity[5], cell_connectivity[6] ],
            [ cell_connectivity[6], cell_connectivity[8] ],
            [ cell_connectivity[8], cell_connectivity[7] ],
            [ cell_connectivity[7], cell_connectivity[5] ],
            [ cell_connectivity[1], cell_connectivity[5] ],
            [ cell_connectivity[2], cell_connectivity[6] ],
            [ cell_connectivity[4], cell_connectivity[8] ],
            [ cell_connectivity[3], cell_connectivity[7] ] ],
            fill(3, 12)
    end
end

function decompose_hexa(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3], cell_connectivity[4] ], 
            [ cell_connectivity[5], cell_connectivity[6], cell_connectivity[7], cell_connectivity[8] ],
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[6], cell_connectivity[5] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[7], cell_connectivity[6] ],
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[8], cell_connectivity[7] ],
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[5], cell_connectivity[8] ] ],
            fill(9, 6)
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[3] ],
            [ cell_connectivity[3], cell_connectivity[4] ], 
            [ cell_connectivity[4], cell_connectivity[1] ],
            [ cell_connectivity[5], cell_connectivity[6] ],
            [ cell_connectivity[6], cell_connectivity[7] ],
            [ cell_connectivity[7], cell_connectivity[8] ],
            [ cell_connectivity[8], cell_connectivity[5] ],
            [ cell_connectivity[1], cell_connectivity[5] ],
            [ cell_connectivity[2], cell_connectivity[6] ],
            [ cell_connectivity[3], cell_connectivity[7] ],
            [ cell_connectivity[4], cell_connectivity[8] ] ], 
            fill(3, 12)
    end
end

function decompose_wedge(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[5], cell_connectivity[4] ], 
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3] ], 
            [ cell_connectivity[2], cell_connectivity[5], cell_connectivity[6], cell_connectivity[3] ],
            [ cell_connectivity[5], cell_connectivity[4], cell_connectivity[6] ],
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[3], cell_connectivity[6] ] ], 
            [9, 5, 9, 5, 9]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[5] ],
            [ cell_connectivity[5], cell_connectivity[4] ], 
            [ cell_connectivity[4], cell_connectivity[1] ],
            [ cell_connectivity[1], cell_connectivity[3] ],
            [ cell_connectivity[2], cell_connectivity[3] ],
            [ cell_connectivity[5], cell_connectivity[6] ],
            [ cell_connectivity[4], cell_connectivity[6] ],
            [ cell_connectivity[3], cell_connectivity[6] ] ], 
            fill(3, 9)
    end
end

function decompose_pyramid(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3], cell_connectivity[4] ], 
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[5] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[5] ], 
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[5] ], 
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[5] ] ], 
            Int[9, 5, 5, 5, 5]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2] ],
            [ cell_connectivity[2], cell_connectivity[3] ],
            [ cell_connectivity[3], cell_connectivity[4] ], 
            [ cell_connectivity[4], cell_connectivity[1] ],
            [ cell_connectivity[1], cell_connectivity[5] ],
            [ cell_connectivity[2], cell_connectivity[5] ],
            [ cell_connectivity[3], cell_connectivity[5] ],
            [ cell_connectivity[4], cell_connectivity[5] ] ], 
            fill(3, 8)
    end
end

function decompose_quadratic_triangle(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[cell_connectivity], Int[22]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[4]], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[5] ], 
            [ cell_connectivity[3], cell_connectivity[1], cell_connectivity[6] ] ], 
            fill(21, 3)
    end
end

function decompose_quadratic_quad(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[cell_connectivity], Int[23]
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[5] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[6] ],
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[7] ],
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[8] ] ],
            fill(21, 4)
    end
end

function decompose_quadratic_tetra(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3], cell_connectivity[5], cell_connectivity[6], cell_connectivity[7] ], 
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[4], cell_connectivity[5], cell_connectivity[9], cell_connectivity[8] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[4], cell_connectivity[6], cell_connectivity[10], cell_connectivity[9] ],
            [ cell_connectivity[3], cell_connectivity[1], cell_connectivity[4], cell_connectivity[7], cell_connectivity[8], cell_connectivity[10] ] ], 
            fill(22, 4)
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[5] ],
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[5] ],
            [ cell_connectivity[3], cell_connectivity[1], cell_connectivity[7] ],
            [ cell_connectivity[1], cell_connectivity[4], cell_connectivity[8] ], 
            [ cell_connectivity[2], cell_connectivity[4], cell_connectivity[8] ], 
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[10] ] ], 
            fill(21, 6)
    end
end

function decompose_quadratic_hexa(cell_connectivity, target::String)
    if target == "Faces"
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[3], cell_connectivity[4], cell_connectivity[9], cell_connectivity[10], cell_connectivity[11], cell_connectivity[12] ],
            [ cell_connectivity[5], cell_connectivity[6], cell_connectivity[7], cell_connectivity[8], cell_connectivity[13], cell_connectivity[14], cell_connectivity[15], cell_connectivity[16] ], 
            [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[6], cell_connectivity[5], cell_connectivity[9], cell_connectivity[18], cell_connectivity[13], cell_connectivity[17] ], 
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[7], cell_connectivity[6], cell_connectivity[10], cell_connectivity[19], cell_connectivity[14], cell_connectivity[18] ],
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[8], cell_connectivity[7], cell_connectivity[11], cell_connectivity[20], cell_connectivity[15], cell_connectivity[19] ],
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[5], cell_connectivity[8], cell_connectivity[12], cell_connectivity[17], cell_connectivity[16], cell_connectivity[20] ] ],
            fill(23, 6)
    else
        return Vector{Int}[ [ cell_connectivity[1], cell_connectivity[2], cell_connectivity[9] ],
            [ cell_connectivity[2], cell_connectivity[3], cell_connectivity[10] ],
            [ cell_connectivity[3], cell_connectivity[4], cell_connectivity[11] ], 
            [ cell_connectivity[4], cell_connectivity[1], cell_connectivity[12] ],
            [ cell_connectivity[5], cell_connectivity[6], cell_connectivity[13] ],
            [ cell_connectivity[6], cell_connectivity[7], cell_connectivity[14] ],
            [ cell_connectivity[7], cell_connectivity[8], cell_connectivity[15] ],
            [ cell_connectivity[8], cell_connectivity[5], cell_connectivity[16] ], 
            [ cell_connectivity[1], cell_connectivity[5], cell_connectivity[17] ],
            [ cell_connectivity[2], cell_connectivity[6], cell_connectivity[18] ],
            [ cell_connectivity[3], cell_connectivity[7], cell_connectivity[19] ],
            [ cell_connectivity[4], cell_connectivity[8], cell_connectivity[20] ] ],
            fill(21, 12)
    end
end

function decompose_cell(dataset::AbstractStaticVTKData, cell_ind::Union{Int, Tuple{Int, Vararg{Int}}}; target::String="")
    decompose_cell(cell_connectivity(dataset, cell_ind), cell_type(dataset, cell_ind))
end
