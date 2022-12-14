module K3surfaces

import Oscar
#import Hecke
import Markdown

using .Oscar
#using .Hecke

import Hecke: @req, orthogonal_submodule, _short_vectors_gram

include("./K3Auto.jl")
include("./Serialize.jl")

function __init__()
  Oscar.registerSerializationType(K3surfaces.Chamber, "Chamber")
  Oscar.registerSerializationType(K3surfaces.BorcherdsCtx, "BorcherdsCtx")
  Hecke.add_assert_scope(:K3Auto)
  Hecke.add_verbose_scope(:K3Auto)
  Hecke.set_verbose_level(:K3Auto, 2)
end

end # module
