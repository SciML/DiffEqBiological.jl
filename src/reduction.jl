using LinearAlgebra

""" 
    canelide(cons_laws::AbstractMatrix, species::AbstractMatrix)

Check if the given species can be simultaneously elided using the given matrix of conservation laws.
"""
function canelide(C::AbstractMatrix, species::AbstractVector) 
    # This amounts to checking that the transformation [ (kept) (elided) ] -> [ (kept) (conserved quantities) ]
    # is injective, equivalently that the linear map (elided) -> (conserved) is injective,
    # equivalently that it has full row rank
    rank(@view C[:,species]) == length(species)  
end

""" 
    getinvmap(C::AbstractMatrix, elided_species::AbstractVector)

Return the inverse transformation from observed species and conserved quantities
to elided species, if it exists (checked via [`canelide`](@ref)), as a pair of matrices
(M_kept, M_cons). The elided species are recovered via 
`state_elided = M_kept * state_kept + M_cons * cons`.
"""
function getinvmap(C::AbstractMatrix, elided_species::AbstractVector)
    @assert canelide(C, elided_species)
    
    n_elided = length(elided_species)
    n_kept = size(C, 2) - n_elided
    n_cons = size(C, 1)
    kept_species = filter(x -> !(x in elided_species), 1:size(C, 2))
    M_cons_float = pinv(C[:,elided_species])
    M_cons = round.(Int, M_cons_float)
    
    @assert M_cons ≈ M_cons_float
    
    M_kept = -M_cons * view(C, :, kept_species)
    (M_kept, M_cons)
end

"""
    conssyms(rn::ReactionSystem, C::AbstractMatrix)

Return a list of new SymbolicUtils symbols for the conserved valued determined by the
conservation law matrix C. These become parameters for the reduced reaction system.

TODO: Use gensymb to ensure uniqueness, or an alternative.
"""
conssyms(rn::ReactionSystem, C::AbstractMatrix) = map(id -> Sym{Number}(id), Symbol("_C", i) for i in 1:size(C,1))

"""
    getsubsdict(rn::ReactionSystem, C::AbstractMatrix, elided_species::AbstractVector, syms_cons))

Return a dictionary of the form `elided_species => f(state_kept, cons)` of symbolic expressions
expressing the elided species in terms of the kept species and the conserved quantities.
"""
function getsubsdict(rn::ReactionSystem, C::AbstractMatrix, elided_species::AbstractVector, syms_cons)
    M_kept, M_cons = getinvmap(C, elided_species)
    ret = Dict{Any,Any}()
    syms = states(rn)
    kept_species_idcs = filter(x -> !(x in elided_species), 1:size(C, 2))
    kept_species = map(idx -> syms[idx], kept_species_idcs)
    for (i, species) in enumerate(elided_species)
        sym = syms[species]
        ret[sym] = sum(view(M_kept, i, :) .* kept_species) + sum(view(M_cons, i, :) .* syms_cons)
    end
    ret
end

"""
    redspecstoich(species, stoich, elided_species)

Removes the elided species from the species and stoichiometry vectors
(can be applied to substrates/products of a reaction).
"""
function redspecstoich(species, stoich, elided_species)
    idcs = findall(x -> !(x in elided_species), species)
    (species[idcs], stoich[idcs])
end

"""
    jumpratelaw(rx; rxvars=get_variables(rx.rate), combinatoric_ratelaw=true)

Variant of `jumpratelaw` that ignores variables not in `rxvars`, for 
partial creation of the full jump rate law. The original doesn't seem to use 
that parameter at all.
"""
function newjumpratelaw(rx; rxvars=get_variables(rx.rate), combinatoric_ratelaw=true)
    ModelingToolkit.@unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = one(eltype(substoich))
        for (i,stoich) in enumerate(substoich)
            substrates[i] in rxvars || continue
            
            s   = substrates[i]
            rl *= s
            isone(stoich) && continue
            for i in one(stoich):(stoich-one(stoich))
                rl *= (s - i)
            end
            combinatoric_ratelaw && (coef *= factorial(stoich))
        end
        
        !isone(coef) && (rl /= coef)
    end
    rl
end

"""
    reducereaction(rx::Reaction, subs::Dict; combinatoric_ratelaw=true)

Creates a reduced version of the reaction using `subs` to eliminate the chosen species.
"""
function reducereaction(rx::Reaction, subs::Dict; combinatoric_ratelaw=true)
    elided_species = keys(subs)
    
    rate_subs = newjumpratelaw(rx; rxvars=elided_species, 
                               combinatoric_ratelaw=combinatoric_ratelaw)
    rate_subs = substitute(rate_subs, subs)
    
    subs_red, substoich_red = redspecstoich(rx.substrates, rx.substoich, elided_species)
    prods_red, prodstoich_red = redspecstoich(rx.products, rx.prodstoich, elided_species)
    netst_red = filter(x -> !(first(x) in elided_species), rx.netstoich)
    
    connection_type = rx.connection_type    
    
    Reaction(rate_subs, subs_red, prods_red, substoich_red, prodstoich_red, netst_red,
             rx.only_use_rate, connection_type)
end

"""
    reducedreactions(rn::ReactionSystem, C::AbstractMatrix, elided_species::AbstractVector, syms_cons; 
                     combinatoric_ratelaw=true)

Creates reduced versions of the reactions in the reaction system without the elided species.
"""
function reducedreactions(rn::ReactionSystem, C::AbstractMatrix, elided_species::AbstractVector, syms_cons; 
                          combinatoric_ratelaw=true)
    subs = getsubsdict(rn, C, elided_species, syms_cons)
    map(rx -> reducereaction(rx, subs; combinatoric_ratelaw=combinatoric_ratelaw), rn.eqs)
end

"""
    reducereacsys(rn::ReactionSystem, elided_species::AbstractVector; combinatoric_ratelaw=true)

Creates reduced reaction system, eliminating the given species simultaneously. Returns a reduced
reaction system and the matrix of conservation laws of the original system.

TODO: Add support for `rn.systems`.
"""
function reducereacsys(rn::ReactionSystem, elided_species::AbstractVector; combinatoric_ratelaw=true)
    S = netstoichmat(rn)
    C = Int.(conservationlaws(S))
    
    @assert isempty(rn.systems)
    
    canelide(C, elided_species) || error("Cannot elide selected species in reaction system")
    
    syms_cons = conssyms(rn, C)
    
    states_new_idcs = filter(x -> !(x in elided_species), 1:size(C, 2))
    states_new = map(idx -> rn.states[idx], states_new_idcs)
    
    ps_new = [ rn.ps; syms_cons ]
    name_new = Symbol(rn.name, "_R")
    
    subs = getsubsdict(rn, C, elided_species, syms_cons)
    observed_new = [ rn.observed; (key ~ value for (key, value) in pairs(subs))... ]

    reactions_new = reducedreactions(rn, C, elided_species, syms_cons; 
                                     combinatoric_ratelaw=combinatoric_ratelaw)
    
    (ReactionSystem(reactions_new, rn.iv, states_new, ps_new, observed_new, 
                    name_new, rn.systems), C)
end
