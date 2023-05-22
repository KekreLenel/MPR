# -------------------------------------------------------------------------
# create_CS_tables.m: produces tables of CS decompositions 
# -------------------------------------------------------------------------
# authors:         Rohan Kekre and Moritz Lenel
# for updates see: https://github.com/KekreLenel/MPR
# -------------------------------------------------------------------------

cs1 = CS_Decomposition[:,ix_bm];
cs0 = CS_Decomposition[:,ix_rnk];

T = ""
T *= @sprintf( "\\%% Real stock return & Data [90\\%% CI] & Model & RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( "Dividend growth news & 33\\%% [-13\\%%,71\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[1]*100, cs0[1]*100);
T *= @sprintf( "\$-\$ Future real rate news & 8\\%% [-6\\%%,21\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[2]*100, cs0[2]*100);
T *= @sprintf( "\$-\$ Future excess return news & 59\\%% [19\\%%,108\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[3]*100, cs0[3]*100);
T *= @sprintf( "\\hline \n");
write(tab_path*"Campbell_Shiller.tex", T);

M1 = M_effects[:,ix_bm];
M0 = M_effects[:,ix_rnk];

T = ""
T *= @sprintf( "  &  Model & RANK & Model/RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( " \$\\Delta\\log(y)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[3], M0[3], M1[3]/M0[3]);
T *= @sprintf( " \$\\Delta\\log(c)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[2], M0[2], M1[2]/M0[2]);
T *= @sprintf( " \$\\Delta\\log(x)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[1], M0[1], M1[1]/M0[1]);
T *= @sprintf( "\\hline \n");
write(tab_path*"m_real_effects.tex", T);

cs1 = CS_Decomposition[:,ix_idio_bm];
cs0 = CS_Decomposition[:,ix_idio_rnk];

T = ""
T *= @sprintf( "\\%% Real stock return & Data [90\\%% CI] & Model & RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( "Dividend growth news & 33\\%% [-13\\%%,71\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[1]*100, cs0[1]*100);
T *= @sprintf( "\$-\$ Future real rate news & 8\\%% [-6\\%%,21\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[2]*100, cs0[2]*100);
T *= @sprintf( "\$-\$ Future excess return news & 59\\%% [19\\%%,108\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[3]*100, cs0[3]*100);
T *= @sprintf( "\\hline \n");
write(tab_path*"Campbell_Shiller_idio.tex", T);

M1 = M_effects[:,ix_idio_bm];
M0 = M_effects[:,ix_idio_rnk];

T = ""
T *= @sprintf( "  &  Model & RANK & Model/RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( " \$\\Delta\\log(y)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[3], M0[3], M1[3]/M0[3]);
T *= @sprintf( " \$\\Delta\\log(c)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[2], M0[2], M1[2]/M0[2]);
T *= @sprintf( " \$\\Delta\\log(x)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f\$ \\\\ \n", M1[1], M0[1], M1[1]/M0[1]);
T *= @sprintf( "\\hline \n");
write(tab_path*"m_real_effects_idio.tex", T);

cs1 = CS_Decomposition[:,ix_interm_bm];
cs0 = CS_Decomposition[:,ix_interm_rnk];


T = ""
T *= @sprintf( "\\%% Real stock return & Data [90\\%% CI] & Model & RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( "Dividend growth news & 33\\%% [-13\\%%,71\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[1]*100, cs0[1]*100);
T *= @sprintf( "\$-\$ Future real rate news & 8\\%% [-6\\%%,21\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[2]*100, cs0[2]*100);
T *= @sprintf( "\$-\$ Future excess return news & 59\\%% [19\\%%,108\\%%] & %6.0f\\%% & %6.0f\\%% \\\\ \n", cs1[3]*100, cs0[3]*100);
T *= @sprintf( "\\hline \n");
write(tab_path*"Campbell_Shiller_interm.tex", T);

M1 = M_effects[:,ix_interm_bm];
M0 = M_effects[:,ix_interm_rnk];

T = ""
T *= @sprintf( "  &  Model & RANK & Model/RANK \\\\ \n");
T *= @sprintf( "\\hline \n");
T *= @sprintf( " \$\\Delta\\log(y)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f \$ \\\\ \n", M1[3], M0[3], M1[3]/M0[3]);
T *= @sprintf( " \$\\Delta\\log(c)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f \$ \\\\ \n", M1[2], M0[2], M1[2]/M0[2]);
T *= @sprintf( " \$\\Delta\\log(x)\$ & \$ %6.0f bp\$ & \$ %6.0f bp\$ & \$ %6.1f \$ \\\\ \n", M1[1], M0[1], M1[1]/M0[1]);
T *= @sprintf( "\\hline \n");
write(tab_path*"m_real_effects_interm.tex", T);

