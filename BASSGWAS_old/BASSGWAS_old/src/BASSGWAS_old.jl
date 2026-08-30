# Copyright (C) 2026 David Helekal
#
# This file is part of BASSGWAS_old.
#
# BASSGWAS_old is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, version 3.
#
# BASSGWAS_old is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with BASSGWAS_old. If not, see <https://www.gnu.org/licenses/>.

module BASSGWAS_old

using LinearAlgebra
using LinearAlgebra.BLAS
using SpecialFunctions
using Distributions
using LogExpFunctions
using LoopVectorization
using Octavian
using Polyester
using Random
using StatsBase
using UnPack

include("utils.jl")
include("ParState.jl")
include("SCache.jl")
include("QCache.jl")
include("likelihoods.jl")
include("bvs_funcs.jl")
export wTGSState2 
include("ProbitSampler.jl")
export ProbitBVSSampler, next_step!, get_state, is_adapted
include("postpred.jl")
export postpred
include("eig_opt.jl")
export optimize_design
end # module BADBwGS
