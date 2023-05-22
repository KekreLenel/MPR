% -------------------------------------------------------------------------
% create_decomp_tables.m: produces decompositions table
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

loc_fileid = fopen([tab_path, filesep, 'decomp_tab_1.tex'],'w');
fprintf(loc_fileid, "& $a$ & $b$ & $c$ \\\\ \\hline\n");
fprintf(loc_fileid, "$d\\log k^i$  &  $%4.0f{bp}$ &  $%4.0f{bp}$ &  $%4.0f{bp}$ \\\\ \\hline\n", 100*100*decomp_mat(:,2,1));
fprintf(loc_fileid, "$n^i/(qk^i)$  &  $%4.1f$ &  $%4.1f$ & $%4.1f$\\\\ \n", decomp_mat(:,3,1));
fprintf(loc_fileid, "$q\\partial (k^i)/\\partial n^i$ &  $%4.1f$ & $%4.1f$ & $%4.1f$\\\\ \n", decomp_mat(:,4,1));
fprintf(loc_fileid, "$d\\log n^i$  &  $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \n", 100*100*decomp_mat(:,5,1));
fprintf(loc_fileid, "$\\partial \\log k^i$  &  $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \\hline \n", 100*100*decomp_mat(:,6,1)+100*100*decomp_mat(:,7,1)); % residual assigned here, bc likely source of error
fclose(loc_fileid);


loc_fileid = fopen([tab_path, filesep, 'decomp_tab_2.tex'],'w');                                                                                                                                                                                            
fprintf(loc_fileid, "& baseline & $\\rho^m = 0.75$ & $\\chi^x = 0$ & $\\chi^w = 0$ \\\\ \\hline\n");                                                                                                                                                        
fprintf(loc_fileid, "$d(\\lambda^a n^a/n)$                                                                          & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \\hline \n", 100*100*decomp_mat2(1, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));   
fprintf(loc_fileid, "$(\\lambda^a/n)\\left(d \\left[w \\ell^a\\right] - (n^a/n)d \\left[w \\ell\\right]\\right)$  & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \n",         100*100*decomp_mat2(2, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));     
fprintf(loc_fileid, "$(\\lambda^a/n)(-(1+i_{-1})/P)(B_{-1}^a  + \\nu^a B_{-1}^g)d \\log p$                                            & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \n",         100*100*decomp_mat2(3, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));     
fprintf(loc_fileid, "$(\\lambda^a/n)(k_{-1}^a - (n^a/n)k_{-1}) d \\pi$                                            & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \n",         100*100*decomp_mat2(4, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));     
fprintf(loc_fileid, "$(\\lambda^a/n)(k_{-1}^a - (n^a/n)k_{-1})(1-\\delta) d q$                                    & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \n",         100*100*decomp_mat2(5, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));     
fprintf(loc_fileid, "$(\\lambda^a/n)d t^a$                                                                        & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ & $%4.0f{bp}$ \\\\ \\hline \n", 100*100*decomp_mat2(6, [ix_bm,ix_rhom,ix_chiX0,ix_chiW0]));     
fclose(loc_fileid); 
