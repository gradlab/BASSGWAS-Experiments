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

#Um :: random effect singular vectors used during sampling
#Sm :: random effect singular values used during sampling
#U :: random effect singular vectors all Ys
#S :: random effect singular values for all Ys
#X :: predictor matrix for all Ys
#idx_obs :: index set of observed Ys
#icept :: intercept draws
#betas :: beta draws
#us :: random effect draws
function postpred(Um, Sm, U, S, X, idx_obs, icept, betas, us, rs; rng=Random.default_rng())

    #P = n_pred
    #N = n_resp
    #M = n_it
    n_it = size(betas, 1)
    n_Y = size(X, 1)

    n_obs = size(idx_obs, 1)
    n_m = size(Sm, 1)
    k = size(S, 1)

    #first extrapolate random effect
    #pseudoinverse for conditional gaussian
    pinv = Um * Diagonal([inv.(Sm); zeros(max(0, n_m - n_obs))])*Um'
    Vm = Um * Diagonal(sqrt.(Sm))
    #observed random effects
    re_obs = Vm * us'
    
    #restricted full svd square root
    Cu = U[idx_obs, :] * Diagonal(sqrt.(S))
    #sanity check that restricted covariances match
    !(isapprox(Cu*Cu', Vm*Vm'; rtol=1e-6, atol=1e-6)) && 
        error("Covariance doesn't match restricted covariance.")
    
    rD = Diagonal(sqrt.(rs))    
    #scaled iid normal draws
    Zs = rand(rng, Normal(), k, n_it) * rD

    #conditional means for all random effects in the svd basis
    cond_means = Cu' * pinv * re_obs
    #conditional covariance for all random effects in the svd basis
    cond_cov = Diagonal(ones(k)) - Cu'*pinv*Cu 
    #square root of cond cov via svd
    sv = svd(cond_cov)
    B = sv.U * Diagonal(sqrt.(sv.S))
    #non-centred random effect coefficients acting in the space spanned by the svd basis 
    re_coeffs = cond_means .+ (B*Zs) 
    
    #full svd square root
    V = U * Diagonal(sqrt.(S))
    #sanity check that predicted random effects for observations match those that were observed, up to precision
    #the condition number appears pretty bad in practice
    !(isapprox(V[idx_obs, :]*re_coeffs, re_obs; rtol=1e-3, atol=1e-3)) && 
        error("Random effects of observations don't match. $(maximum(abs.(V[idx_obs, :]*re_coeffs .- re_obs)))")


    K = Matrix([X ones(n_Y) V]') #model matrix [predictors, intercept, random effect] NxP
    bs = Matrix([betas icept re_coeffs']') #PxM
    linpred = bs'*K #bs' * K' #(K*bs)' #MxN

    lp_success = similar(linpred)
    lp_fail = similar(linpred)

    for i in eachindex(linpred, lp_success, lp_fail)
        lp_success[i] = logccdf(Normal(linpred[i],1.0), 0.0)
        lp_fail[i] = logcdf(Normal(linpred[i],1.0), 0.0)
    end

    return lp_success, lp_fail
end
