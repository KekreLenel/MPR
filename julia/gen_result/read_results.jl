raw_grid_jld = load(data_path*"num_params.jld")
n_I                 = raw_grid_jld["n_I"];
n_states            = raw_grid_jld["n_states"];
n_active_dims       = raw_grid_jld["n_active_dims"];
n_interp            = raw_grid_jld["n_interp"];
n_shocks            = raw_grid_jld["n_shocks"];
n_spread            = raw_grid_jld["n_spread"];
smolyak_d           = raw_grid_jld["smolyak_d"];
n_params            = raw_grid_jld["n_prm"];
n_sim_periods       = raw_grid_jld["n_sim_periods"];
n_irf_periods       = raw_grid_jld["n_irf_periods"];

raw_grid_jld = load(data_path*"grid_locs.jld")
k_grid_mean         = raw_grid_jld["k_grid_mean"];
k_grid_dev          = raw_grid_jld["k_grid_dev"];
s_a_grid_mean   = raw_grid_jld["s_a_grid_mean"];
s_a_grid_dev    = raw_grid_jld["s_a_grid_dev"];
s_c_grid_mean   = raw_grid_jld["s_c_grid_mean"];
s_c_grid_dev    = raw_grid_jld["s_c_grid_dev"];
w_grid_mean         = raw_grid_jld["w_grid_mean"];
w_grid_dev          = raw_grid_jld["w_grid_dev"];
dis_grid_mean       = raw_grid_jld["dis_grid_mean"];
dis_grid_dev        = raw_grid_jld["dis_grid_dev"];
m_grid_mean         = raw_grid_jld["m_grid_mean"];
m_grid_dev          = raw_grid_jld["m_grid_dev"];

# there is no chi0_vec in there now, but I hope it will be stored there

# The following lines are inherited from matlab, not sure it works exactly correctly in julia

for i in eachindex(series_list)
    global sid = series_list[i]
    include("read_series.jl");
end