function triangulate_cell(cell_connectivity::Vector{Int}, cell_type::Int)
    if cell_type == 5 || cell_type == 6 || cell_type == 10
        return decompose_cell(cell_connectivity, cell_type, target = "Faces")[1]
    elseif cell_type == 8
        return Vector{Int}[[cell_connectivity[1], cell_connectivity[2], cell_connectivity[4]], 
            [cell_connectivity[1], cell_connectivity[4], cell_connectivity[3]]]
    elseif cell_type == 9
        return Vector{Int}[[cell_connectivity[1], cell_connectivity[2], cell_connectivity[3]],
            [cell_connectivity[1], cell_connectivity[3], cell_connectivity[4]]]
    elseif cell_type == 11 || cell_type == 12 || cell_type == 13 || cell_type == 14
        tris = Vector{Int}[]
        for (_cell, _cell_type) in zip(decompose_cell(cell_connectivity, cell_type, target = "Faces")...)
            append!(tris, triangulate_cell(_cell, _cell_type))
        end
        return tris
    else
        throw("Cell type $cell_type is not supported.")
    end
end

function triangulate_cell_glmesh(cell_connectivity::Vector{Int}, cell_type::Int)
    if cell_type == 5 || cell_type == 6 || cell_type == 10
        return map(Face{3,UInt32, -1}, decompose_cell(cell_connectivity, cell_type, target = "Faces")[1] .- 1)
    elseif cell_type == 8
        return Face{3,UInt32,-1}[Face{3,UInt32,-1}(cell_connectivity[1]-1, cell_connectivity[2]-1, cell_connectivity[4]-1), 
            Face{3,UInt32,-1}(cell_connectivity[1]-1, cell_connectivity[4]-1, cell_connectivity[3]-1)]
    elseif cell_type == 9
        return Face{3,UInt32,-1}[Face{3,UInt32,-1}(cell_connectivity[1]-1, cell_connectivity[2]-1, cell_connectivity[3]-1), 
            Face{3,UInt32,-1}(cell_connectivity[1]-1, cell_connectivity[3]-1, cell_connectivity[4]-1)]
    elseif cell_type == 11 || cell_type == 12 || cell_type == 13 || cell_type == 14
        tris = Face{3,UInt32,-1}[]
        for (_cell, _cell_type) in zip(decompose_cell(cell_connectivity, cell_type, target = "Faces")...)
            append!(tris, triangulate_cell_glmesh(_cell, _cell_type))
        end
        return tris
    else
        throw("Cell type $cell_type is not supported.")
    end
end
