% -------------------------------------------------------------------------
% read_series.m: reads in individual result series  
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

if strcmp(series_id(1:3),'sim') == 1
    n_periods = n_sim_periods; 
else 
    n_periods  = n_irf_periods;
end

fid = fopen([data_path, strcat(series_id, '_state_series.txt')]);
data_text       = textscan(fid, repmat('%q',1,n_periods*smolyak_d),'Headerlines',0);
data_text       = [data_text{:}];
state_series      = reshape(str2double(data_text),[n_periods, smolyak_d]);
fclose(fid);

fid = fopen([data_path, strcat(series_id, '_vars_series.txt')]);
data_text       = textscan(fid, repmat('%q',1,n_periods*n_interp),'Headerlines',0);
data_text       = [data_text{:}];
other_vars_series      = reshape(str2double(data_text),[n_periods, n_interp]);
fclose(fid);

fid = fopen([data_path, strcat(series_id, '_shock_series.txt')]);
data_text       = textscan(fid, repmat('%q',1,n_periods*n_shocks),'Headerlines',0);
data_text       = [data_text{:}];
shock_series      = reshape(str2double(data_text),[n_periods, n_shocks]);
fclose(fid);

eval(strcat(series_id, '_series(:,:,ccc) = [shock_series, state_series, other_vars_series];'));
