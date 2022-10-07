
################################################################################
#
# Serialization
#
################################################################################

############################################################
# Chamber
import Oscar: SerializerState, DeserializerState

function save_internal(s::SerializerState, D::Chamber)
    return Dict(
        :BorcherdsData => save_type_dispatch(s, D.data),
        :weyl_vector => save_type_dispatch(s, D.weyl_vector),
        :walls => save_type_dispatch(s, D.walls),
        :parent_wall => save_type_dispatch(s, D.parent_wall)
    )
end

function load_internal(s::DeserializerState, ::Type{Chamber}, dict::Dict)
    weyl_vector = load_type_dispatch(s, fmpz_mat, dict[:weyl_vector])
    walls = load_type_dispatch(s, Vector{fmpz_mat}, dict[:walls])
    parent_wall = load_type_dispatch(s, fmpz_mat, dict[:parent_wall])
    data = load_type_dispatch(s, BorcherdsData, dict[:BorcherdsData])
    return Chamber(data, weyl_vector, parent_wall, walls)
end

############################################################
# BorcherdsData
function save_internal(s::SerializerState, D::BorcherdsData)
    return Dict(
        :L => save_type_dispatch(s, D.L),
        :S => save_type_dispatch(s, D.S),
        :compute_OR => save_type_dispatch(s, D.compute_OR), #needs to be worked out
    )
end

function load_internal(s::DeserializerState, ::Type{BorcherdsData}, dict::Dict)
    L = load_type_dispatch(s, ZLat, dict[:L])
    S = load_type_dispatch(s, ZLat, dict[:S])
    compute_OR = load_type_dispatch(s, Bool, dict[:compute_OR])

    return BorcherdsData(L, S, compute_OR)
end
