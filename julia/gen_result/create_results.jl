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
using Dates
using Plots
using LaTeXStrings

results_folder = "../output/tmp/res_1/"
src_folder = "../src/"
param_id = 2;

include(src_folder*"SmolyakApprox.jl")
include(src_folder*"utils.jl")
include(src_folder*"paramsets.jl")
include(src_folder*"init.jl")

new_results = 1;
mpr = init(prm_vec[1]);
n_sim_periods = mpr.n_sim_periods; 
n_irf_periods = mpr.n_irf_periods; 
n_periods = n_sim_periods; 

fig_path  = "../output/figures/";
tab_path  = "../output/tables/";

n_comp = 9;
data_path_1 = "../output/tmp/res_1/";
data_path = "";

# Parameterizations
ix_bm         = 1;
ix_rnk        = 2;
ix_rhom       = 3;
ix_chiX0      = 4;
ix_chiW0      = 5;
ix_idio_bm    = 6;
ix_idio_rnk   = 7;
ix_interm_bm  = 8;
ix_interm_rnk = 9;

state_series        = 0.0
other_vars_series   = 0.0
shock_series        = 0.0
temp_series = 0.0
decomp_mat = zeros(3,7,n_comp)
decomp_mat2 = zeros(6,n_comp)
M_effects = zeros(3,n_comp)
CS_Decomposition = zeros(3,n_comp)
N = 200 # 200
collected_irfs = zeros(N-2,14,4,n_comp)
flag = 0;

series_list = ["sim"  "simdis" "none_irf" "g_irf" "p_irf"  "m_irf"]
series_dict = Dict()
irf_idxs = Dict()

for series_id in series_list
    if series_id[1:3]=="sim" 
        global n_periods = n_sim_periods; 
    else 
        global n_periods  = n_irf_periods;
    end
    series_dict[series_id*"_series"] = zeros(n_periods,mpr.smolyak_d+mpr.n_interp+ mpr.n_shocks,n_comp);
end
sid = ""
ccc = 0

# initialization finished

include("extract_irfs.jl")
global io = ""
for i = 1:n_comp
    global ccc = i
    global param_id = i

    global data_path   = "../output/tmp/res_"*string(ccc)*"/";

    if new_results == 1
        include("read_results.jl");  ## to write
        save(data_path*"data.jld","n_I",n_I)
    else
        data = load(data_path*"data.jld","data")
    end

   # Write parameters into results file
   global io= ""
    global io*= "RUN "* string(ccc)* " "* string(Dates.today())* " "*string(Dates.format(now(), "HH:MM:SS"))*" \n\n"
    global io *= "PARAMETRIZATION\n"
    global io *= "-----------------------\n\n"
        field_list = fieldnames((typeof(prm_vec[param_id])))
            for ppp in eachindex(field_list)
                global io *= string(field_list[ppp])*" "*string(getfield(prm_vec[param_id],field_list[ppp]))*"\n"
    end
    global io *= "-----------------------\n"
    global io *= "avrg p "*string(exp(prm_vec[param_id].disast_p + 0.5*prm_vec[param_id].disast_std^2))*"\n"
    global io *= "std p "*string(sqrt((exp(prm_vec[param_id].disast_std^2)-1)*exp(2*prm_vec[param_id].disast_p +prm_vec[param_id].disast_std^2)))*"\n"
    global io *= "-----------------------\n\n"

    if (ccc == 1)
        mu_bm_vec = [prm_vec[param_id].lmbd_a  1.0-prm_vec[param_id].lmbd_a-prm_vec[param_id].lmbd_c prm_vec[param_id].lmbd_c]';  
    end
    
    ## business cycle moments and tables
    global temp_series = series_dict["sim_series"][:,:,ccc];
    include("extract_series.jl");
    include("calc_moments.jl");
    include("create_moment_tables.jl");
    include("calc_decomposition.jl"); 

    temp_series = series_dict["m_irf_series"][:,:,ccc]; 

    global flag = 1;

    include("extract_series.jl");
    include("calc_CampbellShiller.jl");

    global temp_series = series_dict["none_irf_series"][:,:,ccc];
    global collected_irfs[:,:,1,ccc], ~, ~ = extract_irfs();
    global temp_series = series_dict["m_irf_series"][:,:,ccc];
    global collected_irfs[:,:,2,ccc], irf_idxs, irf_titles = extract_irfs();
    global temp_series = series_dict["p_irf_series"][:,:,ccc];
    global collected_irfs[:,:,3,ccc], ~, ~ = extract_irfs();
    global temp_series = series_dict["g_irf_series"][:,:,ccc];
    global collected_irfs[:,:,4,ccc], ~, ~ = extract_irfs();


    global temp_series = series_dict["simdis_series"][:,:,ccc];
    include("extract_series.jl");


    global io *= "\n\n\nMOMENTS WITH DISASTER\n"

    include("calc_moments.jl");

    write(tab_path*"results_"*string(ccc)*".txt", io);

end



include("create_CS_table.jl");


include("make_figure.jl")
start_t = 1;
len     = 20;

include("plot_irfs.jl");

include("create_decomp_tables.jl");
println("ALL FIGURE GENERATED")