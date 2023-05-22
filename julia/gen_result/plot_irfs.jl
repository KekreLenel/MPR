start_t = 1;
len     = 20;


series_to_plot = [irf_idxs["t1y"]; irf_idxs["Er"]; irf_idxs["Eexc"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
squeeze(a) = dropdims(a, dims = tuple(findall(size(a) .== 1)...))
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],2,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end

legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([1,3],[ix_bm, ix_rnk],[ix_bm, ix_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],layout = (1, 3),size=(900,300))
savefig(p,fig_path*"monetary_return_fig_split1"*".pdf")


series_to_plot = [irf_idxs["exc"]; irf_idxs["sa"]; irf_idxs["infl"]; irf_idxs["q"]; irf_idxs["w"]; irf_idxs["l"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],2,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end
legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([2,3],[ix_bm, ix_rnk],[ix_bm, ix_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],layout = (2,3),size=(900,600))
savefig(p,fig_path*"monetary_return_fig_split2"*".pdf")


series_to_plot = [irf_idxs["x"]; irf_idxs["c"]; irf_idxs["y"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],2,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end
legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([1,3],[ix_bm, ix_rnk],[ix_bm, ix_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],layout = (1,3),size=(900,300))
savefig(p,fig_path*"monetary_return_fig_split3"*".pdf")



series_to_plot = [irf_idxs["t1y"]; irf_idxs["Er"]; irf_idxs["Eexc"]; irf_idxs["exc"]; irf_idxs["sa"]; irf_idxs["y"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],2,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end
legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([2,3],[ix_idio_bm, ix_idio_rnk],[ix_idio_bm, ix_idio_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],layout = (2,3),size=(900,600))
savefig(p,fig_path*"idio_compact_fig_bw"*".pdf")

legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([2,3],[ix_interm_bm, ix_interm_rnk],[ix_interm_bm, ix_interm_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],layout = (2,3),size=(900,600))
savefig(p,fig_path*"interm_compact_fig_bw"*".pdf")


## TFP_compact_fig_bw
series_to_plot = [irf_idxs["i"]; irf_idxs["Er"]; irf_idxs["Eexc"]; irf_idxs["exc"]; irf_idxs["sa"]; irf_idxs["y"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],4,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end
legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([2,3],[ix_bm, ix_rnk],[ix_bm, ix_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],layout = (2,3),size=(900,600))
savefig(p,fig_path*"TFP_compact_fig_bw"*".pdf")


## dis_compact_fig_bw
series_to_plot = [irf_idxs["i"]; irf_idxs["Er"]; irf_idxs["Eexc"]; irf_idxs["exc"]; irf_idxs["sa"]; irf_idxs["y"]];
n_series = size(series_to_plot,1);
y_array = Array{Array{Float64},1}(undef,n_series);
x_array = Array{Array{Float64},1}(undef,n_series);
title_array = Dict();
ylabel_array = Dict();
for iii = 1:n_series
    y_array[iii]      = squeeze(collected_irfs[start_t:len,series_to_plot[iii],3,:] );
    x_array[iii]      = repeat((1:len-start_t+1), 1, n_comp);
    title_array[iii]  = irf_titles[series_to_plot[iii]];
    ylabel_array[iii] = "bp";
end
legend_array = ["Model", "RANK"]; legend_loc1 = 3; legend_loc2 = "ne"; 
ps = make_figure([2,3],[ix_bm, ix_rnk],[ix_bm, ix_rnk],x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
p = plot(ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],layout = (2,3),size=(900,600))
savefig(p,fig_path*"dis_compact_fig_bw"*".pdf")