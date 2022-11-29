
################################################################################
#
# Serialization
#
################################################################################

############################################################
# K3Chamber
import Oscar: SerializerState, DeserializerState

function save_internal(s::SerializerState, D::K3Chamber)
    return Dict(
        :BorcherdsCtx => save_type_dispatch(s, D.data),
        :weyl_vector => save_type_dispatch(s, D.weyl_vector),
        :walls => save_type_dispatch(s, D.walls),
        :parent_wall => save_type_dispatch(s, D.parent_wall)
    )
end

function load_internal(s::DeserializerState, ::Type{K3Chamber}, dict::Dict)
    weyl_vector = load_type_dispatch(s, fmpz_mat, dict[:weyl_vector])
    walls = load_type_dispatch(s, Vector{fmpz_mat}, dict[:walls])
    parent_wall = load_type_dispatch(s, fmpz_mat, dict[:parent_wall])
    data = load_type_dispatch(s, BorcherdsCtx, dict[:BorcherdsCtx])
    return K3Chamber(data, weyl_vector, parent_wall, walls)
end

############################################################
# BorcherdsCtx
function save_internal(s::SerializerState, D::BorcherdsCtx)
    return Dict(
        :L => save_type_dispatch(s, D.L),
        :S => save_type_dispatch(s, D.S),
        :compute_OR => save_type_dispatch(s, D.compute_OR), #needs to be worked out
    )
end

function load_internal(s::DeserializerState, ::Type{BorcherdsCtx}, dict::Dict)
    L = load_type_dispatch(s, ZLat, dict[:L])
    S = load_type_dispatch(s, ZLat, dict[:S])
    compute_OR = load_type_dispatch(s, Bool, dict[:compute_OR])

    return BorcherdsCtx(L, S, compute_OR)
end
