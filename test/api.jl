using Catalyst, DiffEqBase, ModelingToolkit, Test

using ModelingToolkit: value

@parameters t k1 k2
@variables S(t) I(t) R(t)
rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])

specset = Set([value(S) =>1, value(I) => 2, value(R) => 3])
@test issetequal(specset, speciesmap(rs))

pset = Set([value(k1) => 1, value(k2) => 2])
@test issetequal(pset, paramsmap(rs))

rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S,I], [I], [1,1], [2])]
rs2 = ReactionSystem(rxs2, t, [R,I,S], [k2,k1])
@test rs == rs2

rs3 = make_empty_network()
@parameters k3 k4
@variables D(t)
addspecies!(rs3, S)
addspecies!(rs3, D)
addparam!(rs3, k3)
addparam!(rs3, k4)
@test issetequal(species(rs3), [S, D])
@test issetequal(params(rs3), [k3, k4])
addreaction!(rs3, Reaction(k3, [S], [D]))
addreaction!(rs3, Reaction(k4, [S,I], [D]))
merge!(rs, rs3)
addspecies!(rs2, S)
addspecies!(rs2, D)
addparam!(rs2, k3)
addparam!(rs2, k4)
addreaction!(rs2, Reaction(k3, [S], [D]))
addreaction!(rs2, Reaction(k4, [S,I], [D]))
@test rs2 == rs

rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
rs3 = make_empty_network()
addspecies!(rs3, S)
addspecies!(rs3, D)
addparam!(rs3, k3)
addparam!(rs3, k4)
addreaction!(rs3, Reaction(k3, [S], [D]))
addreaction!(rs3, Reaction(k4, [S,I], [D]))
rs4 = merge(rs, rs3)
@test rs2 == rs4

rxs = [Reaction(k1*S, [S,I], [I], [2,3], [2]),
       Reaction(k2*R, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
deps = dependents(rxs[2], rs)
@test isequal(deps, [R,I])
@test isequal(dependents(rxs[1], rs), dependants(rxs[1], rs))
addspecies!(rs, S)
@test numspecies(rs) == 3
addspecies!(rs, S, disablechecks=true)
@test numspecies(rs) == 4
addparam!(rs, k1)
@test numparams(rs) == 2
addparam!(rs, k1, disablechecks=true)
@test numparams(rs) == 3


rnmat = @reaction_network begin
       α, S + 2I --> 2I
       β, 3I --> 2R + S
   end α β

smat = [1 2 0;
        0 3 0]
pmat = [0 2 0;
        1 0 2]
@test all(smat .== substoichmat(rnmat))
@test all(pmat .== prodstoichmat(rnmat))

## testing newly added intermediate complexes reaction networks
rn  = Vector{ReactionSystem}(undef,6)
# mass-action non-catalytic
rn[1] = @reaction_network begin
    k₁, 2A --> B
    k₂, A --> C
    k₃, C --> D
    k₄, B + D --> E
end k₁ k₂ k₃ k₄
Z = [2 0 1 0 0 0 0; 0 1 0 0 0 1 0; 0 0 0 1 0 0 0; 0 0 0 0 1 1 0; 0 0 0 0 0 0 1]
B = [-1 0 0 0; 1 0 0 0; 0 -1 0 0; 0 1 -1 0; 0 0 1 0; 0 0 0 -1; 0 0 0 1]
Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0]
@test all(Z .== complex_stoich_matrix(rn[1]))
@test all(B .== complex_incidence_matrix(rn[1]))
@test all(Δ .== complex_outgoing_matrix(rn[1]))
@test all(Z*B .== netstoichmat(rn[1])')

# mass-action rober
rn[2] = @reaction_network begin
    k₁, A --> B
    k₂, B + B --> C + B
    k₃, B + C --> A + C
end k₁ k₂ k₃
Z = [1 0 0 0 1; 0 1 2 1 0; 0 0 0 1 1]
B = [-1 0 0; 1 0 0; 0 -1 0; 0 1 -1; 0 0 1]
Δ = [-1 0 0; 0 0 0; 0 -1 0; 0 0 -1; 0 0 0]
@test all(Z .== complex_stoich_matrix(rn[2]))
@test all(B .== complex_incidence_matrix(rn[2]))
@test all(Δ .== complex_outgoing_matrix(rn[2]))
@test all(Z*B .== netstoichmat(rn[2])')


#  some rational functions as rates
rn[3] = @reaction_network begin
    k₁ , ∅ --> X₁
    ( k₂/(1 + X₁*X₂ + X₃*X₄ ), k₃/(1 + X₁*X₂ + X₃*X₄ )), 2X₁ + X₂ ↔ 3X₃ + X₄
    k₄, X₄ --> ∅
end k₁ k₂ k₃ k₄
Z = [0 1 2 0 0; 0 0 1 0 0; 0 0 0 3 0; 0 0 0 1 1]
B = [-1 0 0 1; 1 0 0 0; 0 -1 1 0; 0 1 -1 0; 0 0 0 -1]
Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
@test all(Z .== complex_stoich_matrix(rn[3]))
@test all(B .== complex_incidence_matrix(rn[3]))
@test all(Δ .== complex_outgoing_matrix(rn[3]))
@test all(Z*B .== netstoichmat(rn[3])')

# repressilator
rn[4]  = @reaction_network begin
   hillR(P₃,α,K,n), ∅ --> m₁
   hillR(P₁,α,K,n), ∅ --> m₂
   hillR(P₂,α,K,n), ∅ --> m₃
   (δ,γ), m₁ ↔ ∅
   (δ,γ), m₂ ↔ ∅
   (δ,γ), m₃ ↔ ∅
   β, m₁ --> m₁ + P₁
   β, m₂ --> m₂ + P₂
   β, m₃ --> m₃ + P₃
   μ, P₁ --> ∅
   μ, P₂ --> ∅
   μ, P₃ --> ∅
end α K n δ γ β μ;
Z = [0 1 0 0 1 0 0 0 0 0; 0 0 1 0 0 1 0 0 0 0; 0 0 0 1 0 0 1 0 0 0; 0 0 0 0 1 0 0 1 0 0; 0 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 0 1 0 0 1]
B = [-1 -1 -1 1 -1 1 -1 1 -1 0 0 0 1 1 1; 1 0 0 -1 1 0 0 0 0 -1 0 0 0 0 0; 0 1 0 0 0 -1 1 0 0 0 -1 0 0 0 0; 
       0 0 1 0 0 0 0 -1 1 0 0 -1 0 0 0; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
Δ = [-1 -1 -1 0 -1 0 -1 0 -1 0 0 0 0 0 0; 0 0 0 -1 0 0 0 0 0 -1 0 0 0 0 0; 0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0; 
       0 0 0 0 0 0 0 -1 0 0 0 -1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]

@test all(Z .== complex_stoich_matrix(rn[4]))
@test all(B .== complex_incidence_matrix(rn[4]))
@test all(Δ .== complex_outgoing_matrix(rn[4]))
@test all(Z*B .== netstoichmat(rn[4])')

#brusselator
rn[5] = @reaction_network begin
   A, ∅ → X
   1, 2X + Y → 3X
   B, X → Y
   1, X → ∅
end A B
Z = [0 1 2 3 0; 0 0 1 0 1]
B = [-1 0 0 1; 1 0 -1 -1; 0 -1 0 0; 0 1 0 0; 0 0 1 0]
Δ = [-1 0 0 0; 0 0 -1 -1; 0 -1 0 0; 0 0 0 0; 0 0 0 0]
@test all(Z .== complex_stoich_matrix(rn[5]))
@test all(B .== complex_incidence_matrix(rn[5]))
@test all(Δ .== complex_outgoing_matrix(rn[5]))
@test all(Z*B .== netstoichmat(rn[5])')


# some rational functions as rates
rn[6] = @reaction_network begin
    (k₁, k₋₁), X₁ + X₂ <--> X₃ + 2X₄
    (k₂/(1 + X₄*X₅ + X₆*X₇), k₋₂/(1 + X₄*X₅ + X₆*X₇)), 3X₄ + X₅ <--> X₆ + X₇
    (k₃/(1 + X₇ + X₈ + X₉ + X₁₀) , k₋₃/(1 + X₇ + X₈ + X₉ + X₁₀ )), 5X₇ + X₈ <--> X₉ + X₁₀
end k₁ k₋₁ k₂ k₋₂ k₃ k₋₃
Z = [1 0 0 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 2 3 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 1 5 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 1]
B = [-1 1 0 0 0 0; 1 -1 0 0 0 0; 0 0 -1 1 0 0; 0 0 1 -1 0 0; 0 0 0 0 -1 1; 0 0 0 0 1 -1]
Δ = [-1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1]
@test all(Z .== complex_stoich_matrix(rn[6]))
@test all(B .== complex_incidence_matrix(rn[6]))
@test all(Δ .== complex_outgoing_matrix(rn[6]))
@test all(Z*B .== netstoichmat(rn[6])')
