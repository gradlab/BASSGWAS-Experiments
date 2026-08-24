using DataFrames
using CSV
using MCMCChains

cpars_draws = DataFrame(CSV.File("cpars_draws.csv"));
pip_draws = Array{Float64}(DataFrame(CSV.File("pips_draws.csv"))[:,3:end]);
Z = cpars_draws[:,:Z];
Z = Z./sum(Z);
pips = sum(pip_draws .* Z;dims=1)[:];
display(reverse(sort(pips))[1:10])
display(reverse(sortperm(pips))[1:10])
describe(Chains(permutedims(reshape(Array{Float64}(cpars_draws)[:,3:end],1000,:,6),(1,3,2)),[:n, :sigma, :r, :c, :icept, :Z]))


