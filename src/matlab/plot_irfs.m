% -------------------------------------------------------------------------
% plot_irfs.m: creates IRF figures
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

start_t = 1;
len     = 20;


series_to_plot = [irf_idxs.t1y; irf_idxs.Er; irf_idxs.Eexc];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),2,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end
try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([1,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.25 8.5 2.3]);
print(fig_handle,[fig_path, '/monetary_fig_split1'],'-depsc', '-r600','-loose')
fig_handle = make_figure([1,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.25 8.5 2.3]);
print(fig_handle,[fig_path, '/monetary_fig_split1_bw'],'-depsc', '-r600','-loose')
catch me    
end

series_to_plot = [irf_idxs.exc; irf_idxs.sa; irf_idxs.infl; irf_idxs.q; irf_idxs.w; irf_idxs.l];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),2,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end
try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/monetary_fig_split2'],'-depsc', '-r600','-loose')
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/monetary_fig_split2_bw'],'-depsc', '-r600','-loose')
catch me
end

series_to_plot = [irf_idxs.x; irf_idxs.c; irf_idxs.y];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),2,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end
try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([1,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.25 8.5 2.3]);
print(fig_handle,[fig_path, '/monetary_fig_split3'],'-depsc', '-r600','-loose')
fig_handle = make_figure([1,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.25 8.5 2.3]);
print(fig_handle,[fig_path, '/monetary_fig_split3_bw'],'-depsc', '-r600','-loose')
catch me
end


series_to_plot = [irf_idxs.t1y; irf_idxs.Er; irf_idxs.Eexc; irf_idxs.exc; irf_idxs.sa;  irf_idxs.y];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),2,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([2,3],{ix_idio_bm, ix_idio_rnk},{ix_idio_bm, ix_idio_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/idio_compact_fig'],'-depsc', '-r600','-loose')
fig_handle = make_figure([2,3],{ix_idio_bm, ix_idio_rnk},{ix_idio_bm, ix_idio_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/idio_compact_fig_bw'],'-depsc', '-r600','-loose')
catch me 
end
try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([2,3],{ix_interm_bm, ix_interm_rnk},{ix_interm_bm, ix_interm_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/interm_compact_fig'],'-depsc', '-r600','-loose')
fig_handle = make_figure([2,3],{ix_interm_bm, ix_interm_rnk},{ix_interm_bm, ix_interm_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/interm_compact_fig_bw'],'-depsc', '-r600','-loose')
catch me 
end


series_to_plot = [irf_idxs.i; irf_idxs.Er; irf_idxs.Eexc; irf_idxs.exc; irf_idxs.sa;  irf_idxs.y];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),4,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/TFP_compact_fig'],'-depsc', '-r600','-loose')
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/TFP_compact_fig_bw'],'-depsc', '-r600','-loose')
catch me 
end

series_to_plot = [irf_idxs.i; irf_idxs.Er; irf_idxs.Eexc; irf_idxs.exc; irf_idxs.sa;  irf_idxs.y];
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),3,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
legend_array = {"Model", "RANK"}; legend_loc1 = 3; legend_loc2 = "ne"; 
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/dis_compact_fig'],'-depsc', '-r600','-loose')
fig_handle = make_figure([2,3],{ix_bm, ix_rnk},{ix_bm, ix_rnk},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0]);
print(fig_handle,[fig_path, '/dis_compact_fig_bw'],'-depsc', '-r600','-loose')
catch me
end


