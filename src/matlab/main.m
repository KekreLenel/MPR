% -------------------------------------------------------------------------
% main.m: main result file, loops over comparative calibrations and 
% calls functions to read results and to produce tables and figures
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

clear all; close all;

new_results = 1;

% Parameterizations
ix_bm         = 1;
ix_rnk        = 2;
ix_rhom       = 3;
ix_chiX0      = 4;
ix_chiW0      = 5;
ix_idio_bm    = 6;
ix_idio_rnk   = 7;
ix_interm_bm  = 8;
ix_interm_rnk = 9;


addpath ../src/matlab/ 
fig_path  = ['..', filesep, 'output', filesep, 'figures', filesep];
tab_path  = ['..', filesep, 'output', filesep, 'tables'];

fid = fopen(['..',  filesep, 'output', filesep, 'tmp', filesep, 'n_comp.txt'], 'r');
n_comp = fscanf(fid,'%u');

for ccc = 1:n_comp

        data_path   = ['..', filesep, 'output', filesep, 'tmp', filesep, 'res_', num2str(ccc), filesep];
        data_path_1 = ['..', filesep, 'output', filesep, 'tmp', filesep, 'res_1', filesep];
        param_file  = ['..', filesep, 'src', filesep, 'params', filesep, 'param_file_', num2str(ccc), '.csv'];
    
        if new_results == 1

            read_results;
            save([data_path, 'data.mat'])
            
        else

            load([data_path, 'data.mat'])

        end

       %% Write parameters into results file
       fileID = fopen([tab_path, filesep, 'results_', num2str(ccc), '.txt'],'w');
       fprintf(fileID,['RUN ', num2str(ccc), ' ', date, ' ', datestr(now, 'HH:MM:SS'), ' \n\n']);
       fprintf(fileID,'PARAMETRIZATION\n');
       fprintf(fileID,'-----------------------\n\n');
       for ppp = 1:n_params
           fprintf(fileID,'%-16s   %8.4f\n', targ_prm_names{ppp}, targ_prm_table(ppp));
       end
       fprintf(fileID,'-----------------------\n');
       fprintf(fileID,'%-16s   %8.4f\n', 'avrg p', exp(prms.disast_p + 0.5*prms.disast_std^2));
       fprintf(fileID,'%-16s   %8.4f\n', 'std p',  sqrt((exp(prms.disast_std^2)-1)*exp(2*prms.disast_p + prms.disast_std^2)));
       fprintf(fileID,'-----------------------\n\n');

       if (ccc == 1)
           mu_bm_vec = [prms.lmbd_a, 1.0-prms.lmbd_a - prms.lmbd_c, prms.lmbd_c, ]';  
       end
       
       %% business cycle moments and tables
       temp_series = sim_series(:,:,ccc);
       extract_series;

       calc_moments;
       
       create_moment_tables;
       
       calc_decomposition; 

       temp_series = m_irf_series(:,:,ccc); 
       extract_series;
       calc_CampbellShiller;

       [collected_irfs(:,:,1,ccc), ~, ~]                 = extract_irfs(none_irf_series(:,:,ccc), prms);
       [collected_irfs(:,:,2,ccc), irf_idxs, irf_titles] = extract_irfs(m_irf_series(:,:,ccc),    prms);
       [collected_irfs(:,:,3,ccc), ~, ~]                 = extract_irfs(p_irf_series(:,:,ccc),    prms);
       [collected_irfs(:,:,4,ccc), ~, ~]                 = extract_irfs(g_irf_series(:,:,ccc),    prms);

       temp_series = simDis_series(:,:,ccc);
       extract_series;
       fprintf(fileID,'\n\n\nMOMENTS WITH DISASTER\n');
       calc_moments;
       
fclose(fileID);

end

create_CS_tables;
 
plot_irfs;
 
create_decomp_tables;

