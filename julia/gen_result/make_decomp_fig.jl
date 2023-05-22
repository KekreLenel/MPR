using Printf
decomp_mat3 = ones(3,6,2)
decomp_mat4 = ones(6,13)
tab_path  = "../output/tables/";


## table 1
slist = ["\$d\\log k^i \$" 
        "\$ n^i/(qk^i)\$" 
        "\$\\partial (qk^i)/\\partial n^i\$" 
        "\$d\\log n^i\$" 
        "\$\\partial \\log k^i\$"]
T = @sprintf("\\multicolumn{4}{l}{\$- d\\log k/d\\log m = %4.2f \\%% \$} \\\\ \\hline\n",100*decomp_mat3[1,1,1]);
T *= "& \$a\$ & \$b\$ & \$c\$ \\\\ \\hline\n";
for i = 1:length(slist)
    global T *= slist[i]*@sprintf(" &  %4.2f \\%% &  %4.2f \\%% &  %4.2f \\%% \\\\ \\hline\n",decomp_mat3[1,i+1,1],decomp_mat3[3,i+1,1],decomp_mat3[2,i+1,1])
end
write(tab_path*"decomp_tab_1.tex", T);

#3 table 2
slist = ["\$d(n^a/n)\$ " 
        "\$(1/n)\\left(d \\left[w \\ell^a\\right] - (n^a/n)d \\left[w \\ell\\right]\\right)\$ " 
        "\$(1/n)(-(1+i_{-1})/p)b_{-1}^a d \\log p\$" 
        "\$(1/n)(k_{-1}^a - (n^a/n)k_{-1}) d \\pi\$" 
        "\$(1/n)(k_{-1}^a - (n^a/n)k_{-1})(1-\\delta) d q\$" 
        "\$(1/n)d t^a\$"]
T = "& baseline & \$\\rho^m = 0.75\$ & \$\\chi^x = 0\$ & \$\\chi^w = 0\$ \\\\ \\hline\n"
for i = 1:length(slist)
    global T *= slist[i]*@sprintf(" &  %4.2f \\%% &  %4.2f \\%% &  %4.2f \\%%  &  %4.2f \\%% \\\\ \\hline\n",
    100*decomp_mat4[i,1],100*decomp_mat4[i,4],100*decomp_mat4[i,5],100*decomp_mat4[i,6])
end
write(tab_path*"decomp_tab_2.tex", T);

#3 table 3
T = "& baseline & \$\\varphi_w = 0, \\ \\varphi_p=50\$ & \$\\varphi_w = 0, \\ \\varphi_p=50\$ &  \$\\varphi_w = 0, \\ \\varphi_p=0\$\\\\ \\hline\n"
for i = 1:length(slist)
    global T *= slist[i]*@sprintf(" &  %4.2f \\%% &  %4.2f \\%% &  %4.2f \\%%  &  %4.2f \\%% \\\\ \\hline\n",
    100*decomp_mat4[i,1],100*decomp_mat4[i,12],100*decomp_mat4[i,13],100*decomp_mat4[i,6])
end
write(tab_path*"decomp_tab_3.tex", T);


blue        = [0       0.4470   0.7410];
lightblue   = [0       0.9      1.0];
darkblue    = [0       0.2      0.5];
lightred    = [0.9500  0.4250   0.2];
red         = [0.8500  0.3250   0.0980];
yellow      = [0.70    0.640    0.1250];
green       = [0.1     0.75     0.2];
grey        = [0.4     0.4      0.4];

