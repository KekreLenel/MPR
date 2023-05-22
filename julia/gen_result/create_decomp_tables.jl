# -------------------------------------------------------------------------
# create_decomp_tables.m: produces decompositions table
# -------------------------------------------------------------------------
# authors:         Rohan Kekre and Moritz Lenel
# for updates see: https://github.com/KekreLenel/MPR
# -------------------------------------------------------------------------

T = ""
T *= @sprintf("& \$a\$ & \$b\$ & \$c\$ \\\\ \\hline\n");
T *= @sprintf("\$d\\log k^i\$  &  \$%4.0f{bp}\$ &  \$%4.0f{bp}\$ &  \$%4.0f{bp}\$ \\\\ \\hline\n", 100*100*decomp_mat[1,2,1],100*100*decomp_mat[2,2,1],100*100*decomp_mat[3,2,1]);
T *= @sprintf("\$n^i/(qk^i)\$  &  \$%4.1f\$ &  \$%4.1f\$ & \$%4.1f\$\\\\ \n", decomp_mat[1,3,1],decomp_mat[2,3,1],decomp_mat[3,3,1]);
T *= @sprintf("\$q\\partial (k^i)/\\partial n^i\$ &  \$%4.1f\$ & \$%4.1f\$ & \$%4.1f\$\\\\ \n", decomp_mat[1,4,1],decomp_mat[2,4,1],decomp_mat[3,4,1]);
T *= @sprintf("\$d\\log n^i\$  &  \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \n", 100*100*decomp_mat[1,5,1],100*100*decomp_mat[2,5,1],100*100*decomp_mat[3,5,1]);
T *= @sprintf("\$\\partial \\log k^i\$  &  \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \\hline \n", 100*100*decomp_mat[1,6,1]+100*100*decomp_mat[1,7,1],
100*100*decomp_mat[2,6,1]+100*100*decomp_mat[2,7,1],100*100*decomp_mat[3,6,1]+100*100*decomp_mat[3,7,1]); # residual assigned here, bc likely source of error
write(tab_path*"decomp_tab_1.tex", T);

T = ""                                                                                                                                                                                           
T *= @sprintf("& baseline & \$\\rho^m = 0.75\$ & \$\\chi^x = 0\$ & \$\\chi^w = 0\$ \\\\ \\hline\n");                                                                                                                                                        
T *= @sprintf("\$d(\\lambda^a n^a/n)\$                                                                          & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \\hline \n", 100*100*decomp_mat2[1, ix_bm],100*100*decomp_mat2[1, ix_rhom],100*100*decomp_mat2[1, ix_chiX0],100*100*decomp_mat2[1, ix_chiW0]);   
T *= @sprintf("\$(\\lambda^a/n)\\left(d \\left[w \\ell^a\\right] - (n^a/n)d \\left[w \\ell\\right]\\right)\$  & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \n",         100*100*decomp_mat2[2, ix_bm],100*100*decomp_mat2[2, ix_rhom],100*100*decomp_mat2[2, ix_chiX0],100*100*decomp_mat2[2, ix_chiW0]);     
T *= @sprintf("\$(\\lambda^a/n)(-(1+i_{-1})/P)(B_{-1}^a  + \\nu^a B_{-1}^g)d \\log p\$                                            & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \n",         100*100*decomp_mat2[3, ix_bm],100*100*decomp_mat2[3, ix_rhom],100*100*decomp_mat2[3, ix_chiX0],100*100*decomp_mat2[3, ix_chiW0]);     
T *= @sprintf("\$(\\lambda^a/n)(k_{-1}^a - (n^a/n)k_{-1}) d \\pi\$                                            & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \n",         100*100*decomp_mat2[4, ix_bm],100*100*decomp_mat2[4, ix_rhom],100*100*decomp_mat2[4, ix_chiX0],100*100*decomp_mat2[4, ix_chiW0]);     
T *= @sprintf("\$(\\lambda^a/n)(k_{-1}^a - (n^a/n)k_{-1})(1-\\delta) d q\$                                    & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \n",         100*100*decomp_mat2[5, ix_bm],100*100*decomp_mat2[5, ix_rhom],100*100*decomp_mat2[5, ix_chiX0],100*100*decomp_mat2[5, ix_chiW0]);     
T *= @sprintf("\$(\\lambda^a/n)d t^a\$                                                                        & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ & \$%4.0f{bp}\$ \\\\ \\hline \n", 100*100*decomp_mat2[6, ix_bm],100*100*decomp_mat2[6, ix_rhom],100*100*decomp_mat2[6, ix_chiX0],100*100*decomp_mat2[6, ix_chiW0]);     
write(tab_path*"decomp_tab_2.tex", T);