# -------------------------------------------------------------------------
# read_series.jl: reads in individual result series  
# -------------------------------------------------------------------------
# authors:         Rohan Kekre and Moritz Lenel
# for updates see: https://github.com/KekreLenel/MPR
# -------------------------------------------------------------------------


if sid[1:3]=="sim" 
    n_periods = n_sim_periods; 
else 
    n_periods  = n_irf_periods;
end

raw_jld = load(joinpath(data_path, string(sid*"_series", ".jld")))
# change the file path when actually running the code

state_series = raw_jld["state_series"] #dimension is n_period * smolyak_d (7)
other_vars_series = raw_jld["other_vars_series"] # n_period * n_interp
shock_series = raw_jld["shock_series"] # n_period * n_shocks

global series_dict[sid*"_series"][:,:,ccc] = [shock_series state_series other_vars_series];
