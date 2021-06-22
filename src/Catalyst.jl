"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using Reexport, ModelingToolkit
using ModelingToolkit: Symbolic, value, istree, get_states, get_ps, get_iv, get_systems, get_eqs
const DEFAULT_IV = (@parameters t)[1]
@reexport using ModelingToolkit
import MacroTools
import Base: (==), merge!, merge
using Symbolics
using Latexify, Requires
using AbstractAlgebra

# as used in Catlab
const USE_GV_JLL = Ref(false)
function __init__()
    @require Graphviz_jll="3c863552-8265-54e4-a6dc-903eb78fde85" begin
      USE_GV_JLL[] = true
      let cfg = joinpath(Graphviz_jll.artifact_dir, "lib", "graphviz", "config6")
        if !isfile(cfg)
          Graphviz_jll.dot(path -> run(`$path -c`))
        end
      end
    end
  end

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

# reaction network macro
export @reaction_network, @add_reactions

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")

# functions to query network properties
include("networkapi.jl")
export species, params, reactions, speciesmap, paramsmap, numspecies, numreactions, numparams
export make_empty_network, addspecies!, addparam!, addreaction!
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export reaction_complexes, reaction_rates, complex_stoich_matrix, complex_incidence_matrix, complex_outgoing_matrix
export conservationlaws, conservedquantities

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs
import Base.Iterators: flatten
import DataStructures: OrderedDict
import Parameters: @with_kw_noshow
include("graphs.jl")
export Graph, savegraph

end # module
