export weyl_vector, K3Auto, common_invariant, separating_hyperplanes, walls, fingerprint, adjacent_chamber, aut, has_zero_entropy
import Oscar: rays
# set_assert_level(:K3Auto, 0)
# set_verbose_level(:K3Auto, 2)


################################################################################
# Types
################################################################################

mutable struct BorcherdsData
  L::ZLat
  S::ZLat
  SS::ZLat
  R::ZLat
  deltaR::Vector{fmpz_mat}
  dualDeltaR::Vector{fmpz_mat}
  prRdelta::Vector{Tuple{fmpq_mat,fmpq}}
  membership_test
  gramL::fmpz_mat
  gramS::fmpz_mat
  prS::fmpq_mat
  compute_OR::Bool
  # TODO: Store temporary variables for the computations
  # in order to make the core-functions adjacent_chamber and walls
  # as non-allocating as possible.
end

@doc Markdown.doc"""
    BorcherdsData(L::ZLat, S::ZLat, compute_OR::Bool=true) -> BorcherdsData

If `compute_OR` is `false`, then `G` is the subgroup of the orthogonal group of `S`
acting as $\pm 1$ on the discriminant group.
If `compute_OR` is `true`, then `G` consists the subgroup consisting of
isometries of `S` that can be extended to isometries of `L`.
"""
function BorcherdsData(L::ZLat, S::ZLat, compute_OR::Bool=true)
  # transform L to have the standard basis
  # we assume that the basis of L is obtained by completing a basis of R
  # hence we can throw away the R coordinates of a weyl vector when projecting to S
  R = Hecke.orthogonal_submodule(L, S)

  if !iszero(basis_matrix(R)[1:end,1:rank(S)])
    # the following completes the basis of R to a basis of L
    R = lll(Hecke.orthogonal_submodule(L, S))
    basis1 = complete_to_basis(change_base_ring(ZZ,basis_matrix(R)),change_base_ring(ZZ,basis_matrix(L1)))
    basis = vcat(basis1[rank(R)+1:end,:],basis1[1:rank(R),:])
    L2 = lattice(ambient_space(L), basis)

    # Assure that L has the standard basis.
    L3 = Zlattice(gram=gram_matrix(L2))
    V = ambient_space(L3)
    S = lattice(V, basis_matrix(S) * inverse_basis_matrix(L2))

    L = L3
  else
    L1 = Zlattice(gram=gram_matrix(L))
    V = ambient_space(L1)
    S = lattice(V, basis_matrix(S) * inverse_basis_matrix(L))
    L = L1
  end

  SS = Zlattice(gram=gram_matrix(S))

  # precomputations
  R = lll(Hecke.orthogonal_submodule(L, S))
  @assert iszero(basis_matrix(R)[1:end,1:rank(S)])
  bSR = vcat(basis_matrix(S),basis_matrix(R))
  ibSR = inv(bSR)
  I = identity_matrix(QQ,degree(L))
  # prS: L --> S^\vee given with respect to the standard basis of L and the basis of S
  prS = ibSR*I[:,1:rank(S)]#*basis_matrix(S)
  @assert prS[rank(S)+1,:]==0


  if compute_OR
    dd = diagonal(gram_matrix(R))
    @vprint :K3Auto 2 "computing orthogonal group\n"
    OR = orthogonal_group(R)
    @vprint :K3Auto 2 "done\n"
    DR = discriminant_group(R)
    ODR = orthogonal_group(DR)
    imOR = [ODR(hom(DR,DR,[DR(lift(d)*f) for d in gens(DR)])) for f in gens(OR)]
    DS = discriminant_group(S)
    DSS = discriminant_group(SS)
    ODSS = orthogonal_group(DSS)
    orderimOR = order(sub(ODR,imOR)[1])
    @vprint :K3Auto 1 "[O(S):G] = $(order(ODSS)//orderimOR)\n"
    if order(ODR)== orderimOR
      membership_test = (g->true)
    else
      phiSS_S = hom(DSS,DS,[DS(lift(x)*basis_matrix(S)) for x in gens(DSS)])
      phi,i,j = glue_map(L,S,R)
      phi = phiSS_S*inv(i)*phi*j
      img,_ = sub(ODSS,[ODSS(phi*hom(g)*inv(phi)) for g in imOR])
      ds = degree(SS)
      membership_test = (g->ODSS(hom(DSS,DSS,[DSS(vec(matrix(QQ, 1, ds, lift(x))*g)) for x in gens(DSS)])) in img)
    end
  else
    membership_test(g) = is_in_G(SS,g)
  end

  d = exponent(discriminant_group(S))
  Rdual = dual(R)
  sv = short_vectors(rescale(Rdual,-1), 2)
  # not storing the following for efficiency
  # append!(sv,[(-v[1],v[2]) for v in sv])
  # but for convenience we include zero
  T = typeof(sv).parameters[1].parameters[1].parameters[1]
  push!(sv,(zeros(T, rank(Rdual)), QQ(0)))
  rkR = rank(R)
  prRdelta = [(matrix(QQ, 1, rkR, v[1])*basis_matrix(Rdual),v[2]) for v in sv]
  gramL = change_base_ring(ZZ,gram_matrix(L))
  gramS = change_base_ring(ZZ,gram_matrix(S))
  deltaR = [change_base_ring(ZZ,matrix(QQ, 1, rkR, v[1])*basis_matrix(R)) for v in short_vectors(rescale(R,-1),2)]
  dualDeltaR = [gramL*transpose(r) for r in deltaR]
  return BorcherdsData(L, S, SS, R, deltaR, dualDeltaR, prRdelta, membership_test,gramL,gramS,prS,compute_OR)
end




mutable struct Chamber
  weyl_vector::fmpz_mat
  # for v in walls, the corresponding half space is defined by the equation
  # x * gram_matrix(S)*v >= 0, further v is primitive in S (and, in contrast to Shimada, not S^\vee)
  walls::Vector{fmpz_mat}
  parent_wall::fmpz_mat # for the spanning tree
  data::BorcherdsData
  function Chamber()
    return new()
  end
end

export Chamber, BorcherdsData

function Chamber(data::BorcherdsData, weyl_vector::fmpz_mat, parent_wall::fmpz_mat)
  D = Chamber()
  D.weyl_vector = weyl_vector
  D.parent_wall = parent_wall
  D.data = data
  return D
end

function Chamber(data::BorcherdsData, weyl_vector::fmpz_mat, parent_wall::fmpz_mat, walls::Vector{fmpz_mat})
  D = Chamber()
  D.weyl_vector = weyl_vector
  D.parent_wall = parent_wall
  D.data = data
  D.walls = walls
  return D
end

# needed to create sets of Chambers
function Base.hash(C::Chamber)
  return hash(C.weyl_vector[:,1:rank(C.data.S)])
end

# Two chambers are equal if and only if their Weyl vectors
# project to the same point in S
# By the choice of our coordinates this projection is determined
# by the first rank(S) coordinates.
function Base.:(==)(C::Chamber,D::Chamber)
  @req C.data===D.data "must be in the same space"
  return C.weyl_vector[:,1:rank(C.data.S)] == D.weyl_vector[:,1:rank(D.data.S)]
end

@doc Markdown.doc"""
    walls(D::Chamber) -> Vector{fmpz_mat}

Return the walls of the chamber `D`.

The corresponding half space of the wall defined by `v` in `walls(D)` is
$\{x \in S \otimes \RR |  <x,v> \geq 0\}$. Note that `v` is given with respect
to the basis of `S` and is primitive in `S`.
"""
function walls(D::Chamber)
  if !isdefined(D, :walls)
    D.walls = _walls_of_chamber(D.data, D.weyl_vector)
    @assert length(D.walls)>=rank(D.data.S) "$(D.weyl_vector)"
  end
  return D.walls
end

function Oscar.rays(D::Chamber)
  r = reduce(vcat, walls(D), init=zero_matrix(ZZ,0,rank(D.data.SS)))
  rQ = change_base_ring(QQ,r)*gram_matrix(D.data.SS)
  C = positive_hull(rQ)
  Cd = polarize(C)
  L = rays(Cd)
  Lq = fmpq_mat[matrix(QQ,1,rank(D.data.SS),i) for i in L]
  # clear denominators
  Lz = fmpz_mat[change_base_ring(ZZ,i*denominator(i)) for i in Lq]
  # primitive in S
  Lz = fmpz_mat[divexact(i,gcd(vec(i))) for i in Lz]
  @hassert :K3Auto 2 all(all(x>=0 for x in vec(r*gram_matrix(D.data.SS)*transpose(i))) for i in Lz)
  return Lz
end

function Base.show(io::IO, c::Chamber)
  if isdefined(c,:walls)
    print(IOContext(io, :compact => true), "Chamber  in dimension $(length(walls(c)[1])) with $(length(walls(c))) walls")
  else
    print(IOContext(io, :compact => true), "Chamber: $(c.weyl_vector[1,1:rank(c.data.S)])")
  end
end

@doc Markdown.doc"""
    fingerprint(D::Chamber)

Return the hash of the fingerprint of this chamber.

The fingerprint is an invariant computed from the rays and their inner products.
"""
function fingerprint(D::Chamber)
  v = sum(walls(D))
  G = D.data.gramS
  m1 = (v*G*transpose(v))[1,1]
  m2 = [(a*G*transpose(a))[1,1] for a in walls(D)]
  sort!(m2)
  m3 = [(v*G*transpose(a))[1,1] for a in walls(D)]
  sort!(m3)
  m4 = fmpz[]
  for i in 1:length(walls(D))
    for j in 1:i-1
      push!(m4,(walls(D)[i]*G*transpose(walls(D)[j]))[1,1])
    end
  end
  sort!(m4)
  V = Dict{Tuple{fmpz,fmpz},Vector{fmpz_mat}}()
  for w in walls(D)
    i =  (v*G*transpose(w))[1,1]
    j =  (w*G*transpose(w))[1,1]
    if (i,j) in keys(V)
      push!(V[(i,j)],w)
    else
      V[(i,j)] = [w]
    end
  end
  #=
  # so far m5 was not needed to separate the O(S)-orbits
  m5 = []
  for i in keys(V)
    vi = sum(V[i])
    push!(m5, [i,sort!([(vi*G*transpose(j))[1,1] for j in walls(D)])])
  end
  sort!(m5)
  # So far we have only O(S)-invariants. There are also ways to produce G-invariants
  # by working the the images of the rays in the discriminant group and their
  # G-orbits. Perhaps one has to switch to S^\vee primitive vectors in this case.
  =#
  return hash((m1, m2, m3, m4))
end


################################################################################
# close vector functions
################################################################################

@doc Markdown.doc"""
    quadratic_triple -> Vector{Tuple{Vector{Int},fmpq}}

Return $\{x \in Z^n : x Q x^T + 2xb^T + c <=0\}$.

Input:
- `Q` - positive definite matrix
- `b` - vector
- `c` - rational number
"""
function quadratic_triple(Q, b, c; algorithm=:short_vectors, equal=false)
  if algorithm == :short_vectors
    L, p, dist = Hecke._convert_type(Q, b, QQ(c))
    #@vprint :K3Auto 1 ambient_space(L), basis_matrix(L), p, dist
    if equal
      cv = Hecke.close_vectors(L, vec(p), dist, dist, check=false)
    else
      cv = Hecke.close_vectors(L, vec(p), dist, check=false)
    end
  end
  return cv
end

@doc Markdown.doc"""
    short_vectors_affine

Return $\{x \in S : x^2=d, x.v=\alpha \}$.

- `v` - row vector with $v^2 > 0$

Algorithm 2.2 in [Shimada]
"""
function short_vectors_affine(S::ZLat, v::MatrixElem, alpha::fmpq, d)
  gram = gram_matrix(S)
  tmp = v*gram_matrix(ambient_space(S))*transpose(basis_matrix(S))
  v_S = solve_left(gram_matrix(S),tmp)
  sol = short_vectors_affine(gram, v_S, alpha, d)
  B = basis_matrix(S)
  return [s*B for s in sol]
end

function short_vectors_affine(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
  # find a solution <x,v> = alpha with x in L if it exists
  w = gram*transpose(v)
  tmp = FakeFmpqMat(w)
  wn = numerator(tmp)
  wd = denominator(tmp)
  b, x = can_solve_with_solution(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
  if !b
    return fmpq_mat[]
  end
  _, K = left_kernel(wn)
  # (x + y*K)*gram*(x + y*K) = x gram x + 2xGKy + y K G K y

  # now I want to formulate this as a cvp
  # (x +y K) gram (x+yK) ==d
  # (x
  GK = gram*transpose(K)
  Q = K * GK
  b = transpose(x) * GK
  c = (transpose(x)*gram*x)[1,1] - d
  # solve the quadratic triple
  Q = change_base_ring(QQ, Q)
  b = change_base_ring(QQ, transpose(b))
  cv = quadratic_triple(-Q, -b,-QQ(c),equal=true)
  xt = transpose(x)
  cv = [xt+matrix(ZZ,1,nrows(Q),u[1])*K for u in cv]
  @hassert :K3Auto 1 all((v*gram*transpose(u))[1,1]==alpha for u in cv)
  @hassert :K3Auto 1 all((u*gram*transpose(u))[1,1]== d for u in cv)
  return cv #[u for u in cv if (u*gram*transpose(u))[1,1]==d]
end


@doc Markdown.doc"""
    separating_hyperplanes

Return ${x in S : x^2=d, x.v>0, x.h<0}$.

- `S` - a hyperbolic lattice
- `d` - a negative integer
- `v`,`h` - vectors of positive square
"""
function separating_hyperplanes(S::ZLat, v::fmpq_mat, h::fmpq_mat, d)
  V = ambient_space(S)
  @hassert :K3Auto 1 inner_product(V,v,v)[1,1]>0
  @hassert :K3Auto 1 inner_product(V,h,h)[1,1]>0
  gram = gram_matrix(S)
  B = basis_matrix(S)
  vS = solve_left(B,v)
  hS = solve_left(B,h)
  return [a*B for a in separating_hyperplanes(gram,vS,hS,d)]
end

function separating_hyperplanes(gram::fmpq_mat, v::fmpq_mat, h::fmpq_mat, d)
  L = Zlattice(gram=gram)
  n = ncols(gram)
  ch = QQ((h*gram*transpose(h))[1,1])
  cv = QQ((h*gram*transpose(v))[1,1])
  b = basis_matrix(L)
  prW = reduce(vcat,[b[i,:] - (b[i,:]*gram*transpose(h))*ch^-1*h for i in 1:n])
  W = lattice(ambient_space(L), prW, isbasis=false)
  bW = basis_matrix(W)
  # set up the quadratic triple for SW
  gramW = gram_matrix(W)
  s = solve_left(bW, v*prW) * gramW
  Q = gramW + transpose(s)*s*ch*cv^-2

  @vprint :K3Auto 5 Q
  LQ = Zlattice(gram=-Q*denominator(Q))
  S_W = [x[1] for x in short_vectors(LQ,  abs(d*denominator(Q)))]
  append!(S_W,[-x for x in S_W])
  push!(S_W, 0*S_W[1])
  #S_W = quadratic_triple(Q, zero_matrix(QQ,n-1,1), d)
  S_W = [matrix(ZZ,1,nrows(Q),x)*bW for x in S_W]
  S = fmpq_mat[]
  h = change_base_ring(QQ,h)
  for rp in S_W
    rho = abs(d - (rp*gram*transpose(rp))[1,1])*ch^-1
    t,rho = issquare_with_sqrt(rho)
    if !t
      continue
    end
    r = rho*h + rp
    if denominator(r)==1 && (r*gram*transpose(h))[1,1]>0 && (r*gram*transpose(v))[1,1] < 0
      push!(S,r)
    end
  end
  return S
end

@doc Markdown.doc"""
    find_basis(row_matrices::Vector, dim)

Return the first `dim` linearly independent vectors in row_matrices.

We assume that row_matrices consists of row vectors.
"""
function find_basis(row_matrices::Vector, dim::Integer)
  @req length(row_matrices)>=dim > 0 "must contain at least a single vector"
  r = row_matrices[1]
  n = ncols(r)
  B = zero_matrix(base_ring(r), 0, n)
  rk = 0
  for r in row_matrices
    Br = vcat(B, r)
    rk = rank(Br)
    if rk > nrows(B)
      B = Br
    end
    if rk == dim
      break
    end
  end
  @assert rk == dim
  return B
end

find_basis(row_matrices::Vector) = find_basis(row_matrices, ncols(row_matrices[1]))

"""
    is_in_G(S::ZLat, g::fmpz_mat) -> Bool

Return whether the isometry `g` of `S` acts as `+-1` on the discriminant group.
"""
function is_in_G(S::ZLat, g::fmpz_mat)
  D = discriminant_group(S)
  imgs = [D(vec(matrix(QQ,1,rank(S),lift(d))*g)) for d in gens(D)]
  return all(imgs[i] == gens(D)[i] for i in 1:length(gens(D))) || all(imgs[i] == -gens(D)[i] for i in 1:length(gens(D)))
  # OD = orthogonal_group(D)
  # g1 = hom(D,D,[D(lift(d)*g) for d in gens(D)])
  # gg = OD(g1)
  # return isone(gg) || gg == OD(-matrix(one(OD)))
end

Hecke.hom(D::Chamber, E::Chamber) = alg319(gram_matrix(D.data.SS), walls(D), walls(E), D.data.membership_test)

function myaut(D)
  G = -D.data.gramS  # negative since the alg somehow assumes that the lengths are positive
  basis = find_basis(walls(D), ncols(G))
  basis_inv = inv(change_base_ring(QQ,basis))
  gram_basis = basis*G*transpose(basis)
  V = [(ZZ.(vec(i*basis_inv)),(i*G*transpose(i))[1,1]) for i in walls(D)]
  C = Hecke.ZLatAutoCtx([gram_basis])
  fl, Csmall = Hecke.try_init_small(C, true, ZZ(-1), true, V)
  if fl
    return Hecke.auto(Csmall)
  end
end

aut(D::Chamber) = hom(D, D)

# worker for hom and aut
function alg319(gram::MatrixElem, raysD::Vector{fmpz_mat}, raysE::Vector{fmpz_mat}, membership_test)
  n = ncols(gram)
  partial_homs = [zero_matrix(ZZ, 0, n)]
  basis = find_basis(raysD, n)
  gram_basis = basis*gram*transpose(basis)
  # breadth first search
  # Since we expect D and E to be isomorphic,
  # a depth first search with early abort would be more efficient.
  # for now this does not seem to be a bottleneck
  for i in 1:n
    @vprint :K3Auto 3 "level $(i), partial homs $(length(partial_homs)) \n"
    partial_homs_new = fmpz_mat[]
    for img in partial_homs
      extensions = fmpz_mat[]
      k = nrows(img)
      gi = gram*transpose(img)
      for r in raysE
        if (r*gram*transpose(r))[1,1] != gram_basis[k+1,k+1] || (k>0 && r*gi != gram_basis[k+1,1:k])
          continue
        end
        # now r has the correct inner products with what we need
        push!(extensions, vcat(img,r))
      end
      append!(partial_homs_new, extensions)
    end
    partial_homs = partial_homs_new
  end
  basisinv = inv(change_base_ring(QQ, basis))
  homs = fmpz_mat[]
  is_in_hom_D_E(fz) = all(r*fz in raysE for r in raysD)
  vE = sum(raysE) # center of mass of the dual cone
  vD = sum(raysD)
  for f in partial_homs
    f = basisinv*f
    if denominator(f)!=1
      continue
    end
    fz = change_base_ring(ZZ,f)
    if !membership_test(fz)
      continue
    end
    if !(vD*fz == vE)
      continue
    end
    # The center of mass is an interior point
    # Further it uniquely determines the chamber and is compatible with homomorphisms
    # This is basically Remark 3.20
    # -> but this is not true for the center of mass of the dual cone
    push!(homs, fz)
  end
  @hassert :K3Auto 1 all(f*gram*transpose(f)==gram for f in homs)
  return homs
end


@doc Markdown.doc"""
  Compute Delta_w

Output:

Tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S`.

Algorithm 5.8 in [Shi]
"""
# legacy function needed for precomputations
function _alg58(L::ZLat, S::ZLat, R::ZLat, prRdelta, w)
  V = ambient_space(L)
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  Rdual = dual(R)
  Sdual = dual(S)
  rkR = rank(R)
  delta_w = fmpq_mat[]
  for c in n_R
    cm = -c
    for (vr0,vsquare) in prRdelta
      if vsquare != cm
        continue
      end
      a0 = inner_product(V,w,vr0)[1,1]
      if c == 0
        VA = [(vr0,a0)]
      else
        VA = [(vr0,a0),(-vr0,-a0)]
      end
      for (vr,a) in VA
        Sdual_na = short_vectors_affine(Sdual, w, 1 - a, -2 - c)
        for vs in Sdual_na
          vv = vs +  vr
          if myin(vec(vv),L)
            push!(delta_w, vs)
          end
        end
      end
    end
  end
  return delta_w
end

function _alg58(L::ZLat, S::ZLat, R::ZLat, w::MatrixElem)
  Rdual = dual(R)
  sv = short_vectors(rescale(Rdual, -1), 2)
  # not storing the following for efficiency
  # append!(sv,[(-v[1],v[2]) for v in sv])
  # but for convenience we include zero
  T = typeof(sv).parameters[1].parameters[1].parameters[1]
  push!(sv,(zeros(T, rank(Rdual)), QQ(0)))
  rkR = rank(R)
  prRdelta = [(matrix(QQ, 1, rkR, v[1])*basis_matrix(Rdual),v[2]) for v in sv]
  return _alg58(L, S, R, prRdelta, w)
end


function _alg58_short_vector(data::BorcherdsData, w::fmpz_mat)
  L = data.L
  V = ambient_space(L)
  S = data.S
  R = data.R
  wS = w*data.prS
  wSL = wS*basis_matrix(S)
  wL = gram_matrix(L)*transpose(w)
  wSsquare = (wS*data.gramS*transpose(wS))[1,1]
  W = lattice(V, wS*basis_matrix(S))
  N = orthogonal_submodule(S, W)
  # W + N + R < L of finite index
  svp_input = Tuple{fmpq,fmpq_mat,fmpq,Int}[]
  for (rR, rRsq) in data.prRdelta
    if rRsq==2
      continue
    end
    @inbounds rwS = (rR*wL)[1,1]
    alpha = 1 - rwS
    usq = alpha^2*wSsquare^-1 - rRsq
    sq = -2 - usq
    push!(svp_input, (alpha, rR,sq,1))
    alpha = 1 + rwS
    usq = alpha^2*wSsquare^-1 - rRsq
    sq = -2 - usq
    push!(svp_input, (alpha, rR, sq, -1))
  end
  @inbounds bounds = unique!([-i[3] for i in svp_input])
  Ndual = dual(N)
  G = -gram_matrix(Ndual)
  d = denominator(G)
  bounds = [i for i in bounds if divides(d,denominator(i))[1]]
  mi = minimum(bounds)
  ma = maximum(bounds)

  svN = Hecke._short_vectors_gram(Hecke.LatEnumCtx, G,mi,ma, Int64)
  result = fmpq_mat[]
  # treat the special case of the zero vector by copy paste.
  if QQ(0) in bounds
    (rN,sqrN) = (zeros(Int64,rank(Ndual)),0)
    rN1 = zero_matrix(ZZ,1,degree(Ndual))
    found1 = false
    found2 = false
    sqrN = QQ(0)
    for (alpha, rR, sq, si) in svp_input
      if sqrN != sq
        continue
      end
      rr = alpha*wSsquare^-1*wSL + si*rR
      r = rr + rN1
      if !found1 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found1 = true
        push!(result, r*data.prS)
        break
      end
      r = rr - rN1
      if !found2 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found2 = true
        push!(result, r*data.prS)
        break
      end
      if found1 && found2
        break
      end
    end
  end
  for (rN, sqrN) in svN
    if !(sqrN in bounds)
      continue
    end
    rN1 = matrix(ZZ,1,rank(Ndual),rN)*basis_matrix(Ndual)
    found1 = false
    found2 = false
    sqrN = -sqrN
    for (alpha, rR, sq, si) in svp_input
      if sqrN != sq
        continue
      end
      rr = alpha*wSsquare^-1*wSL + si*rR
      r = rr + rN1
      if !found1 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found1 = true
        push!(result, r*data.prS)
        break
      end
      r = rr - rN1
      if !found2 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found2 = true
        push!(result, r*data.prS)
        break
      end
      if found1 && found2
        break
      end
    end
  end
  return result
end

@doc Markdown.doc"""
Compute Delta_w

Output:

Tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S^\vee` given in the basis of S.

Algorithm 5.8 in [Shi]
"""
function _alg58_close_vector(data::BorcherdsData, w::fmpz_mat)
  V = ambient_space(data.L)
  S = data.S
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  @hassert :K3Auto 2 basis_matrix(data.L)==1
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  SSdual = dual(data.SS)
  delta_w = fmpq_mat[]
  wS = w*data.prS
  #wS = solve_left(gram_matrix(S),w*gram_matrix(V)*transpose(basis_matrix(S)))
  Vw = data.gramL*transpose(w)
  # since we do repeated cvp in the same lattice
  # we do the preprocessing here

  # collect the cvp inputs to avoid repeated calculation of the same cvp
  cvp_inputs = Dict{Tuple{fmpq,fmpq},Vector{fmpq_mat}}()
  for c in n_R
    kc = -2-c
    cm = -c
    for (vr,vsquare) in data.prRdelta
      if vsquare != cm
        continue
      end
      a = (vr*Vw)[1,1]  # TODO: could be improved by working in R
      key = (1-a,kc)
      if key in keys(cvp_inputs)
        push!(cvp_inputs[key], vr)
      else
        cvp_inputs[key] = [vr]
      end
      if c != 0
        key1 = (1+a, kc)
        if key1 in keys(cvp_inputs)
          push!(cvp_inputs[key1], -vr)
        else
          cvp_inputs[key1] = [-vr]
        end
      end
    end
  end

  # setup
  gram = gram_matrix(SSdual)
  #B = basis_matrix(S)
  #return [s*B for s in sol]

  # much of this code is copied from
  # short_vectors_affine(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
  # to avoid repeated calculation of the same stuff e.g. K and Q
  # find a solution <x,v> = alpha with x in L if it exists
  ww = transpose(wS)
  tmp = FakeFmpqMat(ww)
  wn = numerator(tmp)
  wd = denominator(tmp)
  _, K = left_kernel(wn)
  K = lll!(K)  # perhaps doing this has no gain?
  # (x + y*K)*gram*(x + y*K) = x gram x + 2xGKy + y K G K y

  # now I want to formulate this as a cvp
  # (x +y K) gram (x+yK) ==d
  # (x
  KG = K*gram
  Q = -KG * transpose(K)
  Qi = inv(Q)
  N = Zlattice(gram=Q,check=false)
  V = ambient_space(N)

  #@show sum(length.(values(cvp_inputs)))
  tmp = zero_matrix(QQ,1,rank(SSdual))


  B = basis_matrix(SSdual)
  KB = K*B
  for (alpha, d) in keys(cvp_inputs)
    can_solve_i, x = can_solve_with_solution(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
    if !can_solve_i
      continue
    end
    # looks like premature optimization ...
    x = change_base_ring(QQ, x)
    b = KG*x
    transpose!(tmp, x)
    # c = (transpose(x)*gram*x)[1,1] - d
    c = (tmp*mul!(x,gram,x))[1,1] - d
    #cv = quadratic_triple(Q,-b,-c,equal=true)
    mul!(b, Qi, b)
    #b = Qi*b
    v = vec(b)
    upperbound = inner_product(V,v,v) + c
    # solve the quadratic triple
    cv = close_vectors(N, v, upperbound, upperbound, check=false)
    mul!(tmp,tmp,B)
    #xtB = transpose(x)*B
    Sdual_na1 = [matrix(ZZ, 1, nrows(Q), u)*KB for (u,_) in cv]
    for v in Sdual_na1
      add!(v,v,tmp)
    end
    Sdual_na2 = [vs*basis_matrix(data.S) for vs in Sdual_na1]
    for i in 1:length(Sdual_na1)
      v = Sdual_na2[i]
      for vr in cvp_inputs[(alpha,d)]
        vv =  v +  vr
        if denominator(vv)==1
          push!(delta_w, Sdual_na1[i])
          break # delta_w is a set, hence we may break
        end
      end
    end
  end

  # the chamber should intersect the boundary only at the QQ-rational points
  @hassert :K3Auto 2 rank(S) == rank(reduce(vcat,[s for s in delta_w]))
  return delta_w
end


@doc Markdown.doc"""
    _walls_of_chamber(data::BorcherdsData, weyl_vector)

Return the walls of the L|S chamber induced by `weyl_vector`.

Corresponds Algorithm 5.11 in [Shi] and calls Polymake.
"""
function _walls_of_chamber(data::BorcherdsData, weyl_vector, alg=:short)
  if alg==:short
    walls1 = _alg58_short_vector(data, weyl_vector)
  elseif alg==:close
    walls1 = _alg58_close_vector(data, weyl_vector)
  end
  if length(walls1)==rank(data.S)
    # shortcut which avoids calling Polymake
    d = rank(data.S)
    walls = Vector{fmpz_mat}(undef,d)
    for i in 1:d
      vs = numerator(FakeFmpqMat(walls1[i]))
      g = gcd(vec(vs))
      if g != 1
        vs = divexact(vs, g)
      end
      walls[i] = vs
    end
    return walls
  end
  i = zero_matrix(QQ, 0, degree(data.SS))
  D = reduce(vcat, (v for v in walls1), init=i)
  P = positive_hull(D)
  r = rays(P)
  d = length(r)
  walls = Vector{fmpz_mat}(undef,d)
  for i in 1:d
    v = matrix(QQ, 1, degree(data.SS), r[i])
    # rescale v to be primitive in S
    vs = numerator(FakeFmpqMat(v))
    g = gcd(vec(vs))
    if g!=1
      vs = divexact(vs, g)
    end
    walls[i] = vs
  end
  return walls
end

function is_S_nondegenerate(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = _alg58(L, S, R, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v*G for v in Delta_w]
  i = zero_matrix(QQ,0,degree(S))
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  return ispointed(P)
end

function inner_point_in_S(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = _alg58(L, S, R, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v*G*transpose(basis_matrix(S)) for v in Delta_w]
  i = zero_matrix(QQ,0,rank(S))
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  Pd = polarize(P)
  @hassert :K3Auto 1 is_pointed(Pd)
  h = sum(rays(Pd))
  h = matrix(QQ,1,rank(S),collect(h))*basis_matrix(S)
  @hassert :K3Auto 1 all(0<inner_product(V,v,h)[1,1] for v in Delta_w)
  return h
end


@doc Markdown.doc"""
    unproject_wall(data::BorcherdsData, vS::fmpz_mat)

Return the (-2)-walls of L containing $v^{perp_S}$ but not $S$.

Based on Algorithm 5.13 in [Shi]
"""
function unproject_wall(data::BorcherdsData, vS::fmpz_mat)
  d = gcd(vec(vS*data.gramS))
  v = QQ(1,d)*(vS*basis_matrix(data.S))  # primitive in Sdual
  vsq = QQ((vS*data.gramS*transpose(vS))[1,1],d^2)

  @hassert :K3Auto 1 vsq>=-2
  rkR = rank(data.R)
  Pv = fmpz_mat[]
  for alpha in 1:Int64(floor(sqrt(Float64(-2//vsq))))
    c = 2 + alpha^2*vsq
    alphav = alpha*v
    for (vr,cc) in data.prRdelta
      # probably we could speed up the for loop and compute
      # and compute the result without a membership test
      # by working modulo ZZ^n directly
      if cc != c
        continue
      end
      vv = alphav + vr
      if denominator(vv)==1
        push!(Pv,change_base_ring(ZZ,vv))
      end
      if c != 0
        vv = alphav - vr
        if denominator(vv)==1
          push!(Pv, change_base_ring(ZZ,vv))
        end
      end
    end
  end
  return Pv
end

@doc Markdown.doc"""
Return return the L|S chamber adjacent to `D` via the wall defined by `v`.
"""
function adjacent_chamber(D::Chamber, v)
  gramL = D.data.gramL
  dualDeltaR = D.data.dualDeltaR
  deltaR = D.data.deltaR
  dimL = ncols(gramL)
  Pv = unproject_wall(D.data, v)
  l = length(Pv)
  @hassert :K3Auto 1 length(Pv) == length(unique(Pv))
  a = 1000000
  @label getu
  a = 2*a
  rep = Array{Tuple{Int,fmpq,Bool}}(undef,l+length(dualDeltaR))
  u = matrix(ZZ, 1, dimL, rand(-a:a, dimL))
  Vw = gramL*transpose(D.weyl_vector)
  Vu = gramL*transpose(u)

  z = zero_matrix(ZZ,1,1)
  for i in 1:length(Pv)
    r = Pv[i]
    mul!(z,r, Vw)
    s = (r*Vu)[1,1]
    divexact!(s,z[1,1])
    if any(rep[j][2]==s for j in 1:i-1)
      @goto getu
    end
    rep[i] = (i, s, false)
  end

  for i in 1:length(deltaR)
    r = deltaR[i]
    mul!(z, r, Vw)
    s = (r*Vu)[1,1]
    divexact!(s,z[1,1])
    if any(rep[j][2]==s for j in 1:l+i-1)
      @goto getu
    end
    rep[l+i] = (i, s, true)
  end
  @hassert :K3Auto 2 length(unique([r[2] for r in rep]))==length(rep)
  sort!(rep, by=x->x[2])
  w = deepcopy(D.weyl_vector)
  tmp = zero_matrix(ZZ,ncols(w),1)
  for (i,s,indualDeltaR) in rep
    if indualDeltaR
      # saves a matrix multiplication
      rdual = dualDeltaR[i]
      r = deltaR[i]
      mul!(z, w, rdual)
      addmul!(w, r, z[1,1])
    else
      r = Pv[i]
      #g = (r*gramL*wt)[1,1]
      #w = w + g*r
      transpose!(tmp,w)
      mul!(tmp,gramL,tmp)
      mul!(z, r, tmp)
      addmul!(w,r, z[1,1])
    end
  end
  # both weyl vectors should lie in the positive cone.
  @assert ((D.weyl_vector)*D.data.gramL*transpose(w))[1,1]>0 "$(D.weyl_vector)    $v\n"
  return Chamber(D.data, w, v)
end


function complete_to_basis(B::fmpz_mat,C::fmpz_mat)
  basis = B
  for j in 1:nrows(C)-nrows(B)
    for i in 1:nrows(C)
      c = C[i,:]
      A = vcat(basis,c)
      h = snf(A)
      if h[1:end,1:nrows(h)]==1
        basis = A
        break
      end
    end
  end
  return basis
end

function K3Auto(S::ZLat, n::Integer; kw...)
  @req n in [10,18,26] "n(=$(n)) must be one of 10,18 or 26"
  L, S, weyl = preprocessingK3Auto(S, n)
  A = K3Auto(L,S,weyl; kw...)
  return collect(A[2]),reduce(append!,values(A[3]),init=Chamber[]), collect(A[4])
end

@doc Markdown.doc"""
Compute the automorphism group of a K3

- `w` - initial Weyl vector
"""
function K3Auto(L::ZLat, S::ZLat, w::fmpq_mat; entropy_abort=false, compute_OR=true, max_nchambers=-1)
  data = BorcherdsData(L, S, compute_OR)
  w = change_base_ring(ZZ,w*inverse_basis_matrix(L))
  # for G-sets
  F = FreeModule(ZZ,rank(S))
  # initialization
  chambers = Dict{UInt64,Vector{Chamber}}()
  explored = Set{Chamber}()
  D = Chamber(data, w, zero_matrix(ZZ, 1, rank(S)))
  waiting_list = [D]

  automorphisms = Set{fmpz_mat}()
  rational_curves = Set{fmpz_mat}()

  # gogo
  ncircles = 0
  ntry = 0
  nchambers = 0
  while length(waiting_list) > 0
    ntry = ntry + 1
    if mod(ntry, 100)==0
      @vprint :K3Auto 2 "largest bucket: $(maximum(length(i) for i in values(chambers))) "
      @vprint :K3Auto 1 "buckets: $(length(chambers)) explored: $(nchambers) unexplored: $(length(waiting_list)) generators: $(length(automorphisms))\n"
    end
    D = popfirst!(waiting_list)
    if D in explored
      continue
    end
    # check G-congruence
    fp = fingerprint(D)  # this is the bottleneck - the computation of the walls.
    if !haskey(chambers,fp)
      chambers[fp] = Chamber[]
    end
    is_explored = false
    for E in chambers[fp]
      @vprint :K3Auto 2 "$(D.weyl_vector)    $(E.weyl_vector)\n"
      gg = hom(D, E)
      if length(gg) > 0
        # enough to add a single homomorphism
        if !(gg[1] in automorphisms)
          push!(automorphisms, gg[1])
          if entropy_abort
            C = lattice(rational_span(S),common_invariant(automorphisms)[2])
            d = diagonal(rational_span(C))
            if 0 > maximum(push!([sign(i) for i in d],-1))
              @vprint :K3Auto 1 "entropy abort \n"
              return data, automorphisms, chambers, rational_curves, false
            end
          end
        end
        is_explored = true
        break
      end
    end
    if is_explored
      continue
    end
    push!(chambers[fp], D)
    push!(explored, D)
    nchambers = nchambers+1
    @vprint :K3Auto 3 "new weyl vector $(D.weyl_vector)\n"

    autD = aut(D)
    autD = [a for a in autD if !isone(a)]
    # we need the orbits of the walls only
    if length(autD) > 0
      for f in autD
        push!(automorphisms, f)
      end
      @vprint :K3Auto 1 "Found a chamber with $(length(autD)) automorphisms\n"
      # compute the orbits
      @vprint :K3Auto 2 "computing orbits"
      Omega = [F(v) for v in walls(D)]
      W = gset(matrix_group(autD),Omega)
      vv = F(D.parent_wall)
      wallsDmodAutD = [representative(w).v for w in orbits(W) if !(vv in w)]
    else
      # the minus shouldnt be necessary ... but who knows?
      wallsDmodAutD = (v for v in walls(D) if !(v==D.parent_wall || -v==D.parent_wall))
    end
    # compute the adjacent chambers to be explored
    for v in wallsDmodAutD
      if -2 == (v*gram_matrix(S)*transpose(v))[1,1]
        # v comes from a rational curve
        push!(rational_curves, v)
        continue
      end
      Dv = adjacent_chamber(D, v)
      push!(waiting_list, Dv)
    end
    if max_nchambers != -1 && ntry > max_nchambers
      return data, automorphisms, chambers, rational_curves, false
    end
  end
  @vprint :K3Auto "$(length(automorphisms)) automorphism group generators\n"
  @vprint :K3Auto "$(nchambers) congruence classes of chambers \n"
  @vprint :K3Auto "$(length(rational_curves)) orbits of rational curves\n"
  return data, automorphisms, chambers, rational_curves, true
end

function dist(V::Hecke.QuadSpace, r::fmpq_mat, h1::fmpq_mat, h2::fmpq_mat)
  if inner_product(V,h1-h2,r)!=0
    return inner_product(V, h1, r)[1,1]//inner_product(V, h1 - h2, r)[1,1]
  else
    return PosInf()
  end
end

function chain_reflect(V::Hecke.QuadSpace, h1, h2, w, separating_walls::Vector{fmpq_mat})
  @hassert :K3Auto 1 inner_product(V,h1,h2)[1,1]>0
  @hassert :K3Auto 1 all(inner_product(V,h1,r)[1,1]>=0 for r in separating_walls)
  @hassert :K3Auto 1 all(inner_product(V,h2,r)[1,1]<=0 for r in separating_walls)
  di(r) = dist(V, r, h1, h2)
  #sort!(separating_walls, by=di)
  separating_walls0 = deepcopy(separating_walls)
  while length(separating_walls0)>0
    _,i = findmax(di(r) for r in separating_walls0)
    r = separating_walls0[i]
    deleteat!(separating_walls0,i)
    if inner_product(V,h2,r)[1,1]>0
      continue
    end
    h2 = h2 + inner_product(V, h2, r)*r
    w = w + inner_product(V, w, r)*r
    separating_walls0 = [r for r in separating_walls0 if inner_product(V,h2,r)[1,1]<0]
    # should be decreasing
    # @vprint :K3Auto 1 length([s for s in separating_walls0 if 0>sign(inner_product(V,h2,s)[1,1])])
  end
  # confirm output .... since I did not yet prove this algorithm .. it looks a bit fishy
  @assert all(inner_product(V,h2,r)[1,1]>=0 for r in separating_walls)
  return h2, w
end

# return QQ(D(weyl)\cap S)
function span_in_S(L, S, weyl)
  R = Hecke.orthogonal_submodule(L, S)
  V = ambient_space(L)
  Delta_w = _alg58(L, S, R, weyl)
  G = gram_matrix(V)
  prSDelta_w = [v*G for v in Delta_w]
  @vprint :K3Auto 2 "Ddual given by $(length(prSDelta_w)) rays\n"
  i = zero_matrix(QQ, 0, degree(S))
  R = Hecke.orthogonal_submodule(L, S)
  Ddual = reduce(vcat, prSDelta_w, init=i)
  Ddual = vcat(Ddual, basis_matrix(R))
  Ddual = vcat(Ddual, -basis_matrix(R))
  Ddual = positive_hull(Ddual)
  D = polarize(Ddual)
  @vprint :K3Auto 3 "calculating rays\n"
  gensN = [matrix(QQ, 1, degree(S), v) for v in vcat(Oscar.rays(D),lineality_space(D))]
  gensN = reduce(vcat, gensN, init=i)
  r = Hecke.rref!(gensN)
  gensN = gensN[1:r,:]
  QQDcapS = lattice(V, gensN)
  return QQDcapS
end

function weyl_vector_non_degenerate(L::ZLat, S::ZLat, u0::fmpq_mat, weyl::fmpq_mat, ample0::fmpq_mat, perturbation_factor=1000)
  V = ambient_space(L)
  ample = ample0
  u = u0


  @vprint :K3Auto 2 "calculating separating hyperplanes\n"
  separating_walls = separating_hyperplanes(L, u, ample, -2)
  @vprint :K3Auto 2 "moving ample class\n"
  u, weyl = chain_reflect(V, ample, u, weyl, separating_walls)

  if is_S_nondegenerate(L,S,weyl)
    return weyl, u, ample
  end
  @vprint :K3Auto 2 "calculating QQDcapS\n"
  QQDcapS = span_in_S(L,S,weyl)

  N = Hecke.orthogonal_submodule(L, QQDcapS)
  N = lll(N)
  @vprint :K3Auto 2 "computing the relevant roots\n"
  @vprint :K3Auto 3 "$(gram_matrix(N))\n"

  sv = short_vectors(N^(-1//1), 2)
  @vprint :K3Auto 2 "done\n"
  relevant_roots = [matrix(QQ,1,rank(N),a[1])*basis_matrix(N) for a in sv]
  T = Hecke.orthogonal_submodule(S, QQDcapS)
  if rank(T)==0
    return weyl,u,u
  end
  @vprint :K3Auto 1 "degeneracy dimension of the chamber $(rank(T))\n"
  @label choose_h
  h = perturbation_factor*ample + matrix(QQ,1,rank(T),rand(-10:10,rank(T)))*basis_matrix(T)
  # roots orthogonal to S do not help. Therefore discard them.
  relevant_roots = [r for r in relevant_roots if inner_product(V,basis_matrix(S),r)!=0]
  if any(inner_product(V,h,r)==0 for r in relevant_roots)
    @goto choose_h
  end
  separating = fmpq_mat[r for r in relevant_roots if sign(inner_product(V, h, r)[1,1])*sign(inner_product(V, u, r)[1,1])<0]
  # fix signs
  for i in 1:length(separating)
    r = separating[i]
    if inner_product(V, u, r)[1,1]>0
      separating[i] = -r
    end
    r = separating[i]
  end
  @hassert :K3Auto 1 all(inner_product(V,h,r)[1,1] > 0 for r in separating)
  @hassert :K3Auto 1 all(inner_product(V,u,r)[1,1] < 0 for r in separating)

  u, weyl = chain_reflect(V, h, u, weyl, separating)
  @hassert :K3Auto 1 all(inner_product(V,u,r)[1,1] > 0 for r in separating)

  @assert is_S_nondegenerate(L,S,weyl)
  return weyl, u, h
end

@doc Markdown.doc"""
Embed `S` into an hyperbolic, even unimodular lattice `L` of rank $n$.

If `n` is `26`, then the orthogonal complement $R = S^\perp$ in `L` has a (-2)-vector.
Or an error is produced (does not enumerate the genus of $R$).
"""
function embed_in_unimodular(S::ZLat, n)
  @vprint :K3Auto 1 "computing embedding in L_$(n) \n"
  r = n - rank(S)
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  G = genus(DR, (0, r))  # genus of R
  R = representative(G)
  R = lll(R)
  R = Zlattice(gram=gram_matrix(R))
  if n==26 && maximum(diagonal(gram_matrix(R)))<-2
    @vprint :K3Auto 2 "checking the embedding"
    @hassert :K3Auto 1 minimum(rescale(R,-1))==2
  end
  SR, iS, iR = orthogonal_sum(S, R)
  V = ambient_space(SR)
  S = lattice(V,basis_matrix(S)*iS.matrix)
  R = lattice(V,basis_matrix(R)*iR.matrix)
  DS,_ = normal_form(discriminant_group(S))
  DR = discriminant_group(R)
  DRn,_ = normal_form(rescale(DR,-1))
  gensDS = [lift(x) for x in gens(DS)]
  imgsDS = [lift(x) for x in gens(DRn)]
  glue = reduce(vcat, [matrix(QQ,1,degree(SR),gensDS[i]+imgsDS[i]) for i in 1:length(gensDS)],init=zero_matrix(QQ,0,degree(S)))
  gensL = vcat(basis_matrix(SR), glue)
  L = lattice(V, gensL, isbasis=false)
  @hassert :K3Auto 1 abs(det(L))==1
  @hassert :K3Auto 1 denominator(gram_matrix(L))==1
  return L, S, iS, R, iR
end


@doc Markdown.doc"""
    weyl_vector(L::ZLat, U0::ZLat)

Return a Weyl vector of `L`.

For `L` of signature (1,25) it uses the 24 holy constructions of the Leech
lattice.

Input:

`L` - an even unimodular lattice of signature (1,9), (1,17) or (1,25)
`U0` - a sublattice of `L` with gram matrix `[0 1; 1 -2]`
"""
# we can do the 24 constructions of the leech lattice
function weyl_vector(L::ZLat, U0::ZLat)
  @vprint :K3Auto 1 "computing an initial Weyl vector \n"
  @assert gram_matrix(U0) == QQ[0 1; 1 -2] "$(gram_matrix(U0))"
  V = ambient_space(L)
  U = U0
  R = Hecke.orthogonal_submodule(L,U)
  R0 = Hecke.orthogonal_submodule(L,U)
  if rank(L)==10
    E8 = R0
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    _, T = isisometric(e8, E8, ambient_representation=false)
    E8 = lattice(V, T * basis_matrix(E8))
    B = vcat(basis_matrix(U), basis_matrix(E8))
    Bdual = inv(gram_matrix(V) * transpose(B))
    # this one does not have ample projection
    weyl = QQ[30 1 1 1 1 1 1 1 1 1] * Bdual
    @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 1240
    return weyl, weyl
  elseif rank(L) == 18
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    e8e8,_,_ = orthogonal_sum(e8, e8)
    while true
      R = Hecke.orthogonal_submodule(L,U)
      @vprint :K3Auto 1 "starting isometry test\n"
      isiso, T = isisometric(e8e8, R, ambient_representation=false)
      @vprint :K3Auto 1 "done\n"
      if isiso
        E8E8 = R
        break
      end
      U = U0
      R = R0

      # compute a 2-neighbor
      v = zero_matrix(QQ,1,rank(R))
      while true
        v = matrix(QQ,1,rank(R),rand(0:1,rank(R)))
        if !iszero(v) && mod(numerator((v*gram_matrix(R)*transpose(v))[1,1]),4)==0
          break
        end
      end
      b = change_base_ring(ZZ, v*gram_matrix(R)*transpose(v)*1//4)
      A = change_base_ring(ZZ, gram_matrix(R)*transpose(v))
      b = change_base_ring(GF(2), b)
      A = change_base_ring(GF(2), A)
      x = lift(solve_left(A, b))
      v = (v + 2*x)*basis_matrix(R)
      @hassert :K3Auto 1 mod(inner_product(V,v,v)[1,1], 8)==0
      u = basis_matrix(U)
      f1 = u[1,:]
      e1 = u[2,:] + u[1,:]
      f2 = -inner_product(V, v, v)*1//4*f1 + 2*e1 + v
      @hassert :K3Auto 1 inner_product(V, f2, f2)==0

      e2 = find_section(L, f2)

      #s = change_base_ring(ZZ, basis_matrix(R)*gram_matrix(V)*transpose(f2))
      #e2 = solve_left(s, matrix(ZZ,1,1,[1]))*basis_matrix(R)
      @hassert :K3Auto 2 inner_product(V, f2, e2)[1,1] == 1
      @hassert :K3Auto 2 inner_product(V,e2,e2)[1,1]==-2
      #e2 = e2 - (inner_product(V,e2,e2)[1,1]*(1//2) + 1)*f2
      u = vcat(f2,e2)
      U = lattice(V,u)
      @hassert :K3Auto 1 gram_matrix(U) == QQ[0 1; 1 -2]
    end
    E8E8 = lattice(V, T * basis_matrix(E8E8))
    B = vcat(basis_matrix(U), basis_matrix(E8E8))
    Bdual = inv(gram_matrix(V) * transpose(B))
    # this one does not have ample projection
    weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
    @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 620
    return weyl, weyl
  elseif rank(L)==26
    while true
        R = lll(Hecke.orthogonal_submodule(L,U))
        m = minimum(rescale(R,-1))
        @vprint :K3Auto 1 "found a lattice of minimum $(m) \n"
        if m==4
          # R is isomorphic to the Leech lattice
          fu = basis_matrix(U)[1,:]
          zu = basis_matrix(U)[2,:]
          u0 = 3*fu+zu
          return fu,u0
        end
        @assert m == 2
        #=
        e8 = rescale(root_lattice(:E,8), -1)
        e8e8,_,_ = orthogonal_sum(e8, e8)
        e8e8e8,_,_ = orthogonal_sum(e8e8, e8)
        # the isometry test seems to be expensive sometimes
        isiso,T = isisometric(e8e8e8, R, ambient_representation=false)
        @vprint :K3Auto 2 root_type(R)[2]
        @vprint :K3Auto 2 "\n"
        if isiso
          break
        end
        U = U0
        R = R0
        =#
        leech,v,h = leech_from_root_lattice(rescale(R,-1))
        # the leech lattice is the h-neighbor of R with respect to v
        # with the attached hyperbolic planes this can be engineered to give an isometry
        @hassert :K3Auto 1 mod(inner_product(V,v,v)[1,1],2*h^2)==0
        u = basis_matrix(U)
        f1 = u[1,:]
        e1 = u[2,:] + u[1,:]
        f2 = -inner_product(V, v, v)*1//(2*h)*f1 + h*e1 + v
        @hassert :K3Auto 1 inner_product(V, f2, f2)==0

        e2 = find_section(L,f2)
        u = vcat(f2,e2)
        U = lattice(V,u)
        @hassert :K3Auto 1 gram_matrix(U) == QQ[0 1; 1 -2]
      end
      E8E8E8 = lattice(V, T * basis_matrix(R))
      B = vcat(basis_matrix(U), basis_matrix(E8E8E8))
      Bdual = inv(gram_matrix(V) * transpose(B))
      # this one does not have ample projection
      weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
      @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 0
      s = find_section(L,weyl)
      u0 = 3*weyl + s
      @hassert :K3Auto 1 inner_product(V, u0,u0)[1,1] == 4
      @hassert :K3Auto 1 inner_product(V, u0,weyl)[1,1] == 1
      return weyl, u0
    end
  error("L must be even, hyperbolic unimodular of rank 10,18,26")
end

function find_section(L::ZLat, f)
  V = ambient_space(L)
  g = [abs(i) for i in vec(inner_product(ambient_space(L),f,basis_matrix(L)))]
  if 1 in g
    i = findfirst(x->x==1,g)
    s = basis_matrix(L)[i,:]
    s = sign(inner_product(ambient_space(L),f,s)[1,1])*s
  else
    # search a smallish section using a cvp
    @hassert :K3Auto 1 inner_product(V,f,f)==0
    A = change_base_ring(ZZ,basis_matrix(L)*gram_matrix(V)*transpose(f))
    ss = solve_left(A,identity_matrix(ZZ,1))
    s = ss*basis_matrix(L)
    k, K = left_kernel(A)
    Kl = Zlattice(gram=K*transpose(K))
    # project ss to K
    sK = solve(change_base_ring(QQ,K*transpose(K)),change_base_ring(QQ,K*transpose(ss)))
    a = QQ(1)
    cv = []
    while true
      cv = Hecke.close_vectors(Kl, vec(-sK), a, check=false)
      if length(cv)>0
        break
      end
      a = a+2
    end
    sK = transpose(sK)
    v0 = 0
    for (v,_) in cv
      v = matrix(QQ, 1, rank(Kl), v)
      v1 = v+sK
      aa = (v1*gram_matrix(Kl)*transpose(v1))[1,1]
      if aa < a
        a = aa
        v0 = v
      end
    end
    s = (v0*K + ss)*basis_matrix(L)
  end

  @vprint :K3Auto 2 "found section of size $(s*transpose(s))\n"
  @hassert :K3Auto 1 inner_product(V,s,f)[1,1]==1
  #return s,f
  s = s - (inner_product(V,s,s)[1,1]//2+1) * f
  @hassert :K3Auto 1 inner_product(V,s,s)[1,1]==-2
  return s
end

@doc Markdown.doc"""
    preprocessingK3Auto(S::ZLat, n::Integer)

Return an embedding of `S` into an even unimodular, hyperbolic lattice L of rank
n as well as an `S`-nondegenerate Weyl vector.
"""
function preprocessingK3Auto(S::ZLat, n::Integer; ample=nothing)
  @req n in [10,18,26] "n must be one of 10, 18 or 26"
  # another example
  S = Zlattice(gram=gram_matrix(S))
  L,S,iS, R,iR = embed_in_unimodular(S::ZLat, n)
  V = ambient_space(L)
  if gram_matrix(S)[1:2,1:2] == QQ[0 1; 1 -2]
    U = lattice(V,basis_matrix(S)[1:2,:])
  else
    # find a hyperbolic plane
    G = gram_matrix(L)
    g1,u1 = lll_gram_indef_with_transform(change_base_ring(ZZ,G))
    # apparently we need to run lll two times to produce a zero
    g2,u2 = lll_gram_indef_with_transform(change_base_ring(ZZ,g1))
    u = u2*u1
    @assert g2 == u*G*transpose(u)
    B = u*basis_matrix(L)
    B = vcat(B[1,:],B[end,:]-B[1,:])
    if inner_product(V,B,B) != QQ[0 1; 1 -2]
      # find an isotropic vector
      if gram_matrix(S)[1,1]==0
        v = basis_matrix(S)[1,:]
      else
        b,v = Hecke._isisotropic_with_vector(gram_matrix(L))
        @assert b
        v = matrix(QQ,1,degree(L),v)
        v = v*basis_matrix(L)
      end
      s = find_section(L,v)
      B = vcat(v,s)
    end
    U = lattice(V, B)
  end
  @assert gram_matrix(U) == QQ[0 1; 1 -2]
  weyl, u0 = weyl_vector(L, U)
  #project the weyl vector to S.
  v = u0*gram_matrix(ambient_space(L))*transpose(basis_matrix(S))
  vv = solve_left(gram_matrix(S),v)
  c = (vv*gram_matrix(S)*transpose(vv))[1,1]
  cc = QQ(Int64(floor(sqrt(Float64(c)))))
  vv = 10//cc * vv
  vv = matrix(QQ,1, ncols(vv),[round(i) for i in vv])
  @assert (vv*gram_matrix(S)*transpose(vv))[1,1]>=0
  ntry = 0
  if ample isa Nothing
    #find a random ample vector ... or use a perturbation of the weyl vector?
    @vprint :K3Auto 1 "searching a random ample vector in S\n"
    while true
      ntry = ntry+1
      fudge = floor(ntry//100)
      h = (fudge*vv+matrix(ZZ,1,rank(S),rand(-10-fudge*10:10+fudge*10, rank(S))))*basis_matrix(S)
      # confirm that h is in the interior of a weyl chamber,
      # i.e. check that Q does not contain any -2 vector and h^2>0
      if inner_product(V,h,h)[1,1]<=0
        continue
      end
      @hassert :K3Auto 1 0 < inner_product(V,h,h)[1,1]
      Q = Hecke.orthogonal_submodule(S, lattice(V, h))
      if length(short_vectors(rescale(Q, -1), 2)) == 0
        break
      end
    end
    if inner_product(V,weyl,h)[1,1]<0
      h = -h
    end
  else
    h = ample*basis_matrix(S)
    Q = orthogonal_submodule(S, lattice(V,h))
    @assert length(short_vectors(rescale(Q,-1),2))==0
  end
  weyl1,u,hh = weyl_vector_non_degenerate(L,S,u0, weyl,h)
  return L,S,weyl1#L,S,u0, weyl,weyl1, h
end


function parse_zero_entropy(filename="/home/simon/Dropbox/Math/MyPapers/zero entropy/CandidatesZeroEntropy_elliptic")
  io = open(filename,"r")
  s = read(io, String)
  s = split(s,"\n")
  s = [a for a in s if length(a) > 0 && a[1:1]!="/"]
  s = [split(a," ") for a in s]
  s = [[ZZ(parse(Int64,i)) for i in a] for a in s]
  res = []
  u = ZZ[0 1; 1 -2]
  for g in s
    n = Int64(sqrt(ZZ(length(g))))
    m = matrix(ZZ,n,n,g)
    m = block_diagonal_matrix([u,-m])
    push!(res,m)
  end
  return res
end


################################################################################
# the 23 holy constructions of the leech lattice
################################################################################

function coxeter_number(ADE::Symbol, n)
  if ADE == :A
    return n+1
  elseif ADE == :D
    return 2*(n-1)
  elseif ADE == :E && n == 6
    return 12
  elseif ADE == :E && n == 7
    return 18
  elseif ADE == :E && n == 8
    return 30
  end
end

function highest_root(ADE::Symbol, n)
  if ADE == :A
    w = [1 for i in 1:n]
  elseif ADE == :D
    w = vcat([1,1],[2 for i in 3:n-1])
    w = vcat(w,[1])
  elseif ADE == :E && n == 6
    w = [1,2,3,2,1,2]
  elseif ADE == :E && n == 7
    w = [2,3,4,3,2,1,2]
  elseif ADE == :E && n == 8
    w = [2,4,6,5,4,3,2,3]
  end
  w = matrix(ZZ, 1, n, w)
  g = gram_matrix(root_lattice(ADE,n))
  @hassert :K3Auto 2 all(0<=i for i in collect(w*g))
  @hassert :K3Auto 2 (w*g*transpose(w))[1,1]==2
  return w
end

function _weyl_vector(R::ZLat)
  weyl = matrix(ZZ,1,rank(R),ones(1,rank(R)))*inv(gram_matrix(R))
  return weyl*basis_matrix(R)
end

@doc Markdown.doc"""
    leech_from_root_lattice(niemeier_lattice::ZLat) -> Zlat,fmpq_mat, Integer

Construct the Leech lattice from `niemeier_lattice` by using one of the 23-holy constructions.

Returns a triple `leech_lattice, v, h` where `leech_lattice` is constructed as
`h`-neighbor of `niemeier_lattice` with respect to `v`. `h` is the Coxeter number
of the Niemer lattice.
"""
function leech_from_root_lattice(niemeier_lattice::ZLat)
  # construct the leech lattice from one of the 23 holy constructions in SPLAG
  # we follow Ebeling
  # there seem to be some signs wrong in Ebeling?
  V = ambient_space(niemeier_lattice)
  ADE, ade, RR = root_lattice_recognition_fundamental(niemeier_lattice)
  length(ade)>0 || error("not a niemeier lattice")
  F = basis_matrix(ADE)
  for i in 1:length(ade)
    F = vcat(F, -highest_root(ade[i]...)*basis_matrix(RR[i]))
  end
  rho = sum(_weyl_vector(r) for r in RR)
  h = coxeter_number(ade[1]...)
  @hassert :K3Auto 1 inner_product(V,rho,rho)== 2*h*(h+1)
  @hassert :K3Auto 1 all(h==coxeter_number(i...) for i in ade)
  rhoB = solve_left(basis_matrix(niemeier_lattice),rho)
  v = QQ(1,h)*transpose(rhoB)
  A = Zlattice(gram=gram_matrix(niemeier_lattice))
  c = QQ(2*(1+1//h))
  sv = [matrix(QQ,1,24,vec(v)-i)*basis_matrix(niemeier_lattice) for (i,_) in Hecke.close_vectors(A, vec(v) ,c,c, check=false)]
  @hassert :K3Auto 1 all(inner_product(V,i,i)==2*(1+1//h) for i in sv)
  @hassert :K3Auto 1 length(sv)^2 == abs(det(ADE))
  G = reduce(vcat,(i for i in sv))
  FG = vcat(F,G)
  K = transpose(kernel(matrix(ZZ,ones(Int,1,nrows(FG))))[2])
  B = change_base_ring(QQ,K)*FG
  B = hnf(FakeFmpqMat(B))
  B = QQ(1,B.den)*change_base_ring(QQ,B.num[end-23:end,:])
  @hassert :K3Auto 1 rank(B)==24
  leech_lattice = lattice(V,B)
  @hassert :K3Auto 1 denominator(gram_matrix(leech_lattice))==1
  leech_lattice = lll(leech_lattice)
  @hassert :K3Auto 1 det(leech_lattice)==1
  @hassert :K3Auto 1 minimum(leech_lattice)==4

  T = torsion_quadratic_module(leech_lattice,intersect(leech_lattice,niemeier_lattice))
  @assert length(gens(T))==1 "I just expect this ... but did not really prove it"
  w = transpose(matrix(lift(gens(T)[1])))

  vniemeier = matrix(ZZ, hcat(ones(Int,1,nrows(F)), zeros(Int,1,nrows(G))))
  vleech = matrix(ZZ, hcat(zeros(Int,1,nrows(F)), ones(Int,1,nrows(G))))
  K = transpose(kernel(vcat(vniemeier, vleech))[2])
  return leech_lattice, h*w, h
end


################################################################################
# entropy
################################################################################

function common_invariant(Gamma)
  return left_kernel(reduce(hcat,[g-1 for g in Gamma]))
end

function has_zero_entropy(S; rank_unimod=26)
  #=
  L,S,iS,R,iR = embed_in_unimodular(S,rank_unimod)
  @assert length(short_vectors(rescale(R,-1),2))>0
  V = ambient_space(L)
  U = lattice(V,basis_matrix(S)[1:2, :])
  @hassert :K3Auto 1 det(U)==-1
  weyl,u0 = weyl_vector(L, U)
  #v = matrix(QQ,ones(Int,1,rank(S)-2))*inv(gram_matrix(S)[3:end,3:end])
  #v = denominator(v)*v
  #h = hcat(QQ[0 0 ], v) *basis_matrix(S)  #an ample vector
  u = basis_matrix(U)
  h = zero_matrix(QQ,1,rank(S))
  v = 3*u[1,:] + u[2,:]
  fudge = 1
  nt = 0
  while true
    h = matrix(QQ,1,rank(S)-2,rand(-5:5,rank(S)-2))
    h = hcat(zero_matrix(QQ,1,2),h)*basis_matrix(S)
    b = inner_product(V,h,h)[1,1]
    bb = ZZ(ceil(sqrt(Float64(abs(b)))/2))+fudge
    h = h + bb*v
    @hassert :K3Auto 1 inner_product(V,h,h)[1,1]>0
    # confirm that h is in the interior of a weyl chamber,
    # i.e. check that Q does not contain any -2 vector and h^2>0
    Q = rescale(Hecke.orthogonal_submodule(S, lattice(V, h)),-1)
    Q = lll(Q)
    @vprint :K3Auto 1 "testing ampleness $(inner_product(V,h,h)[1,1])\n"
    sv = short_vectors(Q,2)
    nt = nt+1
    if nt >10
      fudge = fudge+1
      nt = 0
    end
    if length(sv)>0
      @vprint :K3Auto 1 "not ample\n"
      continue
    end
    @vprint :K3Auto 1 "found ample class $(h)\n"
    @vprint :K3Auto 1 "computing an S-non-degenerate weyl vector\n"
    weyl1,u1 = weyl_vector_non_degenerate(L,S,u0,weyl,h)
    if is_S_nondegenerate(L,S,weyl1)
      weyl = weyl1
      break
    end
  end
  @assert is_S_nondegenerate(L,S,weyl)

  @vprint :K3Auto 1 "preprocessing completed \n"
  if preprocessing_only
    return L,S,weyl
  end
  =#
  L, S, weyl = preprocessingK3Auto(S, rank_unimod)
  @vprint :K3Auto 1 "Weyl vector: $(weyl)\n"
  data, K3Autgrp, chambers, rational_curves, _ = K3Auto(L,S,weyl, entropy_abort=false)
  C = lattice(rational_span(S),common_invariant(K3Autgrp)[2])
  d = diagonal(rational_span(C))

  return maximum(push!([sign(i) for i in d],-1)), data, K3Autgrp, chambers, rational_curves
end


function check_zero_entropy(candidate::ZLat, filename="")
  z, data, K3Autgrp, chambers, rational_curves = has_zero_entropy(candidate)
  chambers = reduce(append!,values(chambers),init=Chamber[])
  chambers = [c.weyl_vector for c in chambers]
  save("$(filename).data", [data.L, data.S, collect(K3Autgrp), chambers, collect(rational_curves)])
  io = open(filename, "w")
  println(io, gram_matrix(candidate))
  if z > 0
    println(io, "elliptic")
  elseif z == 0
    println(io, "parabolic")
  elseif z < 0
    println(io, "hyperbolic")
  end
  close(io)
end

function check_zero_entropy(candidates::Vector,postfix="",wa="a")
  ioelliptic = open("elliptic$(postfix)", wa)
  ioparabolic = open("parabolic$(postfix)", wa)
  iohyperbolic = open("hyperbolic$(postfix)", wa)
  close(ioelliptic)
  close(ioparabolic)
  close(iohyperbolic)
  for S in candidates
    e = has_zero_entropy(S)[1]
    ioelliptic = open("elliptic$(postfix)", "a")
    ioparabolic = open("parabolic$(postfix)", "a")
    iohyperbolic = open("hyperbolic$(postfix)", "a")
    if e>0
      println(ioelliptic, gram_matrix(S))
    elseif e==0
      println(ioparabolic, gram_matrix(S))
    elseif e < 0
      println(iohyperbolic, gram_matrix(S))
    end
    close(ioelliptic)
    close(ioparabolic)
    close(iohyperbolic)
  end
end

@attr function inverse_basis_matrix(L::ZLat)
  @hassert :K3Auto 1 degree(L) == rank(L)
  return inv(basis_matrix(L))::fmpq_mat
end

function myin(v, L::ZLat)
  return all(denominator(i)==1 for i in v*inverse_basis_matrix(L))
end
