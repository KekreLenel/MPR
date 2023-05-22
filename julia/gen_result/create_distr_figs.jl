fontsize = 14;
Interpreter = "latex";

blue        = [0       0.4470   0.7410];
lightblue   = [0       0.9      1.0];
darkblue    = [0       0.2      0.5];
lightred    = [0.9500  0.4250   0.2];
red         = [0.8500  0.3250   0.0980];
yellow      = [0.70    0.640    0.1250];
green       = [0.1     0.75     0.2];
grey        = [0.4     0.4      0.4];

h1 = histogram(tht_a_series,nbins=20,colors = blue,alpha = 0.5,xlabel=L"s^a",legend = false,fontsize = fontsize,size=(400,300))
savefig(h1,fig_path*"/sa_distribution_fig"*".pdf")

ps = Dict()
for iii = 1:3
    if iii !=3
        ps[iii] = plot(((iii-1)*400+1):iii*400,tht_a_series[(iii-1)*400+1+5000:iii*400+5000], color = "black",xlims=((iii-1)*400+1, iii*400),ylims=(0.15,0.25),ylabel = L"s^a", linewidth=1,fontsize = fontsize);      
    else
        ps[iii] = plot(((iii-1)*400+1):iii*400,tht_a_series[(iii-1)*400+1+5000:iii*400+5000],  color = "black",xlims=((iii-1)*400+1, iii*400),ylims =(0.15,0.25),ylabel = L"s^a",xlabel=L"quarters", linewidth=1,fontsize = fontsize);
    end
end
p = plot(ps[1],ps[2],ps[3], layout = (3, 1),legend= false,fontsize = fontsize,size=(400,800));
savefig(p,fig_path*"/share_path_fig"*".pdf")
