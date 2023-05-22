% -------------------------------------------------------------------------
% read_results.m: reads in result files produced by fortran code 
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

%% Read in grid specifications and numerical parameters

fid             = fopen([data_path, 'num_params.csv'], 'r');
data_text       = textscan(fid, repmat('%q',1,10),'Headerlines',0);
data_text       = [data_text{:}];
gridsize_vec    = str2double(data_text);
fclose(fid);

fid             = fopen([data_path, 'grid_locs.csv'], 'r');
data_text       = textscan(fid, repmat('%q',1,12),'Headerlines',0);
data_text       = [data_text{:}];
gridlocs_vec    = str2double(data_text);
fclose(fid);

fid             = fopen([data_path, 'extra_data.csv']);
data_text       = textscan(fid, repmat('%q',1,3),'Headerlines',0);
data_text       = [data_text{:}];
extra_data_vec  = str2double(data_text);
chi0_vec        = extra_data_vec(1:3);
fclose(fid);

% Read numerical parameters
n_I                 = gridsize_vec(1);
n_states            = gridsize_vec(2);
n_active_dims       = gridsize_vec(3);
n_interp            = gridsize_vec(4);
n_shocks            = gridsize_vec(5);
n_spread            = gridsize_vec(6);
smolyak_d           = gridsize_vec(7);
n_params            = gridsize_vec(8);
n_sim_periods       = gridsize_vec(9);
n_irf_periods       = gridsize_vec(10);

k_grid_mean         = gridlocs_vec(1);
k_grid_dev          = gridlocs_vec(2);
s_a_grid_mean       = gridlocs_vec(3);
s_a_grid_dev        = gridlocs_vec(4);
s_c_grid_mean       = gridlocs_vec(5);
s_c_grid_dev        = gridlocs_vec(6);
w_grid_mean         = gridlocs_vec(7);
w_grid_dev          = gridlocs_vec(8);
dis_grid_mean       = gridlocs_vec(9);
dis_grid_dev        = gridlocs_vec(10);
m_grid_mean         = gridlocs_vec(11);
m_grid_dev          = gridlocs_vec(12);


%% Read in parameter file

fid = fopen(param_file, 'r');
fmt = repmat('%q ', 1, n_params );
headline = textscan(fid, fmt, 1, 'Delimiter',',','HeaderLines',0);
targ_prm_names = [headline{1:end}];
targ_prm_table = csvread(param_file,1);
n_targ_prm = size(targ_prm_table,2);

prms = cell2struct(cell(1,n_targ_prm),targ_prm_names,2);
for ppp = 1:n_targ_prm
    prms = setfield(prms,targ_prm_names{ppp},targ_prm_table(ppp));
end
prms.debt_to_equity = 0.5;
prms.lmbd_c = 1 - prms.lmbd_a - prms.lmbd_b;


%% Read simulated series and impulse responses

series_id = 'sim';
read_series;

series_id = 'simDis';
read_series;

series_id = 'none_irf';
read_series;

series_id = 'g_irf';
read_series;

series_id = 'p_irf';
read_series;

series_id = 'm_irf';
read_series;




