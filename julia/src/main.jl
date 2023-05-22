# include("pkg_setup.jl")

using Random, Distributions
using Parameters
using FastGaussQuadrature
using SmolyakApprox
using LinearAlgebra
using Profile
using Revise
using UnPack, StructArrays
using Interpolations
using ForwardDiff
using ProgressMeter
using NLsolve
using JLD
using Printf
using BenchmarkTools
using DelimitedFiles

results_folder = "../output/tmp/res_1/"
param_id = 2;

include("SmolyakApprox.jl")
include("utils.jl")
include("paramsets.jl")
include("init.jl")
include("calc_portfolio_share.jl")
include("calc_excess_bond_nom.jl")
include("calc_equilibrium_and_update.jl")
include("calc_bond.jl")
include("calc_val.jl")
include("calc_unexpected_transition.jl")
include("mod_result.jl")
include("mod_decomp.jl")


for i = 1:9 
    println(("="^10)*" paramset $i "*("="^10))
    global param_id = i;
    global results_folder = "../output/tmp/res_"*string(param_id)*"/";
    mpr_vars = init(prm_vec[param_id])
    calc_val(mpr_vars,0,0,i)
    create_results(mpr_vars,i)
    monetary_decomposition(mpr_vars,i)
end