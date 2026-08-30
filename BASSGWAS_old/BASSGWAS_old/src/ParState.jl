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

mutable struct ParState2{Tx, Ty, Tv, Tn, Tr}
    Xaug::Tx
    Y::Ty
    U::Tx
    S::Tv
    p::Tn
    n::Tn
    alpha::Tr
    beta::Tr
    c::Tr
    r::Tr
    s::Tr
end

function ParState2(X, Y, U, S, alpha, beta, c, r, s)
    alpha, beta, c, r, s = promote(alpha, beta, c, r, s)
    p = size(X, 2)    
    n = size(Y, 1)

    Xaug = [X ones(n)]
    #Xaug = [X ones(n)]
    return ParState2(Xaug, Y, U, S, p, n, alpha, beta, c, r, s)
end
