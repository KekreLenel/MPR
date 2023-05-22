function make_figure( fig_size, series_set, series_set_2, x_vec_array, y_vec_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2)

    fontsize = 14;
    
    n_series  = size(series_set,1);
    n_series2 = size(series_set_2,1);
    
    blue        = RGB(0 ,  0.4470 , 0.7410);
    lightblue   = RGB(0     ,  0.9    ,  1.0);
    darkblue    = RGB(0   ,    0.2   ,   0.5);
    lightred    = RGB(0.9500 , 0.4250   ,0.2);
    red         = RGB(0.8500,  0.3250  , 0.0980);
    yellow      = RGB(0.70 ,   0.640 ,   0.1250);
    green       = RGB(0.1  ,   0.75  ,  0.2);
    grey        = RGB(0.4  ,   0.4    ,  0.4);
    
    
    if n_series2 == 1
        color_list = [darkblue];
        style_list = ["-"];
    elseif n_series2 == 2
        color_list = [blue, red];
        style_list = ["-", "-"];
    elseif n_series2 == 3
        color_list = [blue, red, green];
        style_list = ["-", "-","-"];
    elseif n_series2 == 4
        color_list = [blue, red, yellow, green];
        style_list = ["-", "-","-", "-"];
    else
        color_list = [lightblue, blue, darkblue];
        style_list = ["-", "-","-"];
    end
    
    n_subplots = fig_size[1]*fig_size[2];
    
    ps = Dict();
    for nnn = 1:n_subplots 
    
        x_vec = x_vec_array[nnn];
        y_vec = y_vec_array[nnn];
        
        length = size(x_vec,1);
    
        h_cur = hline([0], color =:black, linestyle = :solid, label = false, alpha=0.5, lw=1)
        if (!isempty(legend_array)) && (nnn == legend_loc1)
            for ppp = 1:n_series
                plot!(h_cur,x_vec[:,series_set[ppp]],y_vec[:,series_set[ppp]] ,title=latexstring(title_array[nnn]),ylabel = latexstring(ylabel_array[nnn]),fontsize = fontsize,legendfontsize = fontsize-4,linewidth=2, color= color_list[ppp],Linestyle=:solid,label = latexstring(legend_array[ppp])); # location
            end
        else
            for ppp = 1:n_series
                plot!(h_cur,x_vec[:,series_set[ppp]],y_vec[:,series_set[ppp]] ,title=latexstring(title_array[nnn]),ylabel = latexstring(ylabel_array[nnn]),fontsize = fontsize,legendfontsize = fontsize-4,linewidth=2, color= color_list[ppp],Linestyle=:solid,label = false);
            end
        end
        
        ps[nnn] = h_cur;
        xlims!((x_vec[1,1],x_vec[end,1]));
        if minimum(y_vec[:,series_set_2[:]]) - maximum(y_vec[:,series_set_2[:]]) ==0
            ylims!((-0.01 + maximum(y_vec[:,series_set_2[:]]), maximum(y_vec[:,series_set_2[:]]) + 0.01));    
        elseif isnan(minimum(y_vec[:,series_set_2[:]]))
            ylims!((-0.01 , 0.01))
        else
            ylims!((minimum(y_vec[:,series_set_2[:]]) - 0.05*maximum(abs.(y_vec[:,series_set_2[:]])), maximum(y_vec[:,series_set_2[:]])+ 0.05*maximum(abs.(y_vec[:,series_set_2[:]]))));
        end
    end   
    return ps
end
    
    
    