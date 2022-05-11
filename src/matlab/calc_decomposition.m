% -------------------------------------------------------------------------
% calc_decomposition.m: calculate decomposition results
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

%% 1. Read Decomposition
fid = fopen([data_path_1, 'irf_factor.txt']);
data_text       = textscan(fid, repmat('%q',1,1),'Headerlines',0);
data_text       = [data_text{:}];
irf_factor_1    = reshape(str2double(data_text),[1, 1]);
fclose(fid);

fid = fopen([data_path, 'irf_factor.txt']);
data_text       = textscan(fid, repmat('%q',1,1),'Headerlines',0);
data_text       = [data_text{:}];
irf_factor      = reshape(str2double(data_text),[1, 1]);
fclose(fid);
% -------------------------------------------------- %
% Load State S 0
% -------------------------------------------------- %
fid = fopen([data_path, 'Policies_S_0.txt']);
data_text       = textscan(fid, repmat('%q',1,n_interp),'Headerlines',0);
data_text       = [data_text{:}];
Policy_S_0      = reshape(str2double(data_text),[1, n_interp]);
fclose(fid);

fid = fopen([data_path, 'States_S_0.txt']);
data_text       = textscan(fid, repmat('%q',1,smolyak_d),'Headerlines',0);
data_text       = [data_text{:}];
State_S_0      = reshape(str2double(data_text),[1, smolyak_d]);
fclose(fid);

S_0 = [0,0,0,State_S_0,Policy_S_0];

% -------------------------------------------------- %
% Load State S Star
% -------------------------------------------------- %
fid = fopen([data_path, 'Policies_S_Star.txt']);
data_text       = textscan(fid, repmat('%q',1,n_interp),'Headerlines',0);
data_text       = [data_text{:}];
Policy_S_Star      = reshape(str2double(data_text),[1, n_interp]);
fclose(fid);

fid = fopen([data_path, 'States_S_Star.txt']);
data_text       = textscan(fid, repmat('%q',1,smolyak_d),'Headerlines',0);
data_text       = [data_text{:}];
State_S_Star      = reshape(str2double(data_text),[1, smolyak_d]);
fclose(fid);

S_Star = [0,0,0,State_S_Star,Policy_S_Star];

% -------------------------------------------------- %
% Load State S1
% -------------------------------------------------- %
fid = fopen([data_path, 'Policies_S_1.txt']);
data_text       = textscan(fid, repmat('%q',1,n_interp),'Headerlines',0);
data_text       = [data_text{:}];
Policy_S_1      = reshape(str2double(data_text),[1, n_interp]);
fclose(fid);

fid = fopen([data_path, 'States_S_1.txt']);
data_text       = textscan(fid, repmat('%q',1,smolyak_d),'Headerlines',0);
data_text       = [data_text{:}];
State_S_1      = reshape(str2double(data_text),[1, smolyak_d]);
fclose(fid);

S_1 = [0,0,0,State_S_1,Policy_S_1];

% -------------------------------------------------- %
% Counterfactual Policies
% -------------------------------------------------- %
fid = fopen([data_path, 'Counterfactual_Policies_A.txt']);
data_text       = textscan(fid, repmat('%q',1,3*4),'Headerlines',0);
data_text       = [data_text{:}];
CF_Policies_A      = reshape(str2double(data_text),[3, 4]);
fclose(fid);

fid = fopen([data_path, 'Counterfactual_Policies_B.txt']);
data_text       = textscan(fid, repmat('%q',1,3*4),'Headerlines',0);
data_text       = [data_text{:}];
CF_Policies_B      = reshape(str2double(data_text),[3, 4]);
fclose(fid);

% -------------------------------------------------- %
%% 2. Implement Decomposition
% -------------------------------------------------- %
% (1) Derive the individual policies at S Star
% -------------------------------------------------- %
sav_a_0 = S_Star(jx_sav_a); 
k_a_0   = S_Star(jx_k_a); 
b_a_0   = S_Star(jx_sav_a) - S_Star(jx_k_a)*S_Star(jx_q); 
c_a_0   = prms.lmbd_a * S_Star(jx_c_a);

sav_b_0 = S_Star(jx_sav_b); 
k_b_0   = S_Star(jx_k_b); 
b_b_0   = S_Star(jx_sav_b) - S_Star(jx_k_b)*S_Star(jx_q); 
c_b_0   = prms.lmbd_b * S_Star(jx_c_b);

sav_c_0 = S_Star(jx_sav_c); 
k_c_0   = S_Star(jx_k_c); 
b_c_0   = S_Star(jx_sav_c) - S_Star(jx_k_c)*S_Star(jx_q); 
c_c_0   = prms.lmbd_c * S_Star(jx_c_c);



% -------------------------------------------------- %
% (2) Derive the individual policies at S1
% -------------------------------------------------- %
sav_a_1 = S_1(jx_sav_a); 
k_a_1   = S_1(jx_k_a); 
b_a_1   = S_1(jx_sav_a) - S_1(jx_k_a)*S_1(jx_q); 
c_a_1   = prms.lmbd_a * S_1(jx_c_a);

sav_b_1 = S_1(jx_sav_b); 
k_b_1   = S_1(jx_k_b); 
b_b_1   = S_1(jx_sav_b) - S_1(jx_k_b)*S_1(jx_q); 
c_b_1   = prms.lmbd_b * S_1(jx_c_b);

sav_c_1 = S_1(jx_sav_c); 
k_c_1   = S_1(jx_k_c); 
b_c_1   = S_1(jx_sav_c) - S_1(jx_k_c)*S_1(jx_q); 
c_c_1   = prms.lmbd_c * S_1(jx_c_c);

% -------------------------------------------------- %
% (3) Scale Counterfactual Policies
% -------------------------------------------------- %
c_a_CFA = CF_Policies_A(1,1) * prms.lmbd_a;
c_b_CFA = CF_Policies_A(2,1) * prms.lmbd_b;
c_c_CFA = CF_Policies_A(3,1) * prms.lmbd_c;

b_a_CFA = CF_Policies_A(1,3) * prms.lmbd_a;
b_b_CFA = CF_Policies_A(2,3) * prms.lmbd_b;
b_c_CFA = CF_Policies_A(3,3) * prms.lmbd_c;

k_a_CFA = CF_Policies_A(1,4) * prms.lmbd_a;
k_b_CFA = CF_Policies_A(2,4) * prms.lmbd_b;
k_c_CFA = CF_Policies_A(3,4) * prms.lmbd_c;

c_a_CFB = CF_Policies_B(1,1) * prms.lmbd_a;
c_b_CFB = CF_Policies_B(2,1) * prms.lmbd_b;
c_c_CFB = CF_Policies_B(3,1) * prms.lmbd_c;

b_a_CFB = CF_Policies_B(1,3) * prms.lmbd_a;
b_b_CFB = CF_Policies_B(2,3) * prms.lmbd_b;
b_c_CFB = CF_Policies_B(3,3) * prms.lmbd_c;

k_a_CFB = CF_Policies_B(1,4) * prms.lmbd_a;
k_b_CFB = CF_Policies_B(2,4) * prms.lmbd_b;
k_c_CFB = CF_Policies_B(3,4) * prms.lmbd_c;

% -------------------------------------------------- %
% (4) Compute Direct Effects
% 3 by 3 Table, a, b, c
% -------------------------------------------------- %
Direct_Effects_k = [k_a_CFA-k_a_0, k_b_CFA-k_b_0, k_c_CFA-k_c_0];
Direct_Effects_b = [b_a_CFA-b_a_0, b_b_CFA-b_b_0, b_c_CFA-b_c_0];
Direct_Effects_c = [c_a_CFA-c_a_0, c_b_CFA-c_b_0, c_c_CFA-c_c_0];

Return_Effects_k = [k_a_CFA-k_a_CFB, k_b_CFA-k_b_CFB, k_c_CFA-k_c_CFB];
Return_Effects_b = [b_a_CFA-b_a_CFB, b_b_CFA-b_b_CFB, b_c_CFA-b_c_CFB];
Return_Effects_c = [c_a_CFA-c_a_CFB, c_b_CFA-c_b_CFB, c_c_CFA-c_c_CFB];

% -------------------------------------------------- %
% (5) Compute Total Effects
% -------------------------------------------------- %
Total_Effects_k = [k_a_1-k_a_0, k_b_1-k_b_0, k_c_1-k_c_0];
Total_Effects_b = [b_a_1-b_a_0, b_b_1-b_b_0, b_c_1-b_c_0];
Total_Effects_c = [c_a_1-c_a_0, c_b_1-c_b_0, c_c_1-c_c_0];

% -------------------------------------------------- %
% (6) Load MPK/MPB/MPC at S_Star
% -------------------------------------------------- %
MPK_a = S_Star(jx_MPK1) / S_Star(jx_q);
MPK_b = S_Star(jx_MPK2) / S_Star(jx_q);
MPK_c = S_Star(jx_MPK3) / S_Star(jx_q);
MPC_a = S_Star(jx_MPC1);
MPC_b = S_Star(jx_MPC2);
MPC_c = S_Star(jx_MPC3);
MPB_a = 1 - S_Star(jx_MPC1) - S_Star(jx_MPK1);
MPB_b = 1 - S_Star(jx_MPC2) - S_Star(jx_MPK2);
MPB_c = 1 - S_Star(jx_MPC3) - S_Star(jx_MPK3);

MP_k = [MPK_a, MPK_b, MPK_c];
MP_b = [MPB_a, MPB_b, MPB_c];
MP_c = [MPC_a, MPC_b, MPC_c];

% -------------------------------------------------- %
% (7) Load Net Worth Changes
% -------------------------------------------------- %
n_a_0 = S_Star(jx_n_a); 
n_b_0 = S_Star(jx_n_b); 
n_c_0 = S_Star(jx_n_c); 

aggr_wealth_1 = (S_0(jx_sav_a) + S_0(jx_sav_b) + S_0(jx_sav_c))*(S_1(jx_pi) + (1-prms.delta)*S_1(jx_q));
n_a_1 = prms.s_bar_a*prms.xi*aggr_wealth_1 + (1-prms.xi)*S_0(jx_sav_a)* ... 
((1-S_0(jx_sh_a))*(S_1(jx_pi) + (1-prms.delta)*S_1(jx_q))/S_0(jx_q) + S_0(jx_sh_a)*S_0(jx_nom_i)/S_1(jx_infl));
n_b_1 = (1- prms.s_bar_c - prms.s_bar_a)*prms.xi*aggr_wealth_1 + (1-prms.xi)*S_0(jx_sav_b)* ... 
((1-S_0(jx_sh_b))*(S_1(jx_pi) + (1-prms.delta)*S_1(jx_q))/S_0(jx_q) + S_0(jx_sh_b)*S_0(jx_nom_i)/S_1(jx_infl)); 
n_c_1 = prms.s_bar_c*prms.xi*aggr_wealth_1 + (1-prms.xi)*S_0(jx_sav_c)* ... 
((1-S_0(jx_sh_c))*(S_1(jx_pi) + (1-prms.delta)*S_1(jx_q))/S_0(jx_q) + S_0(jx_sh_c)*S_0(jx_nom_i)/S_1(jx_infl));
% -------------------------------------------------- %
% (8) Indirect Effects and Residuals
% -------------------------------------------------- %
Indirect_Effects_k = [MPK_a*(n_a_1-n_a_0), ...
                      MPK_b*(n_b_1-n_b_0), ...
                      MPK_c*(n_c_1-n_c_0)];
Indirect_Effects_b = [MPB_a*(n_a_1-n_a_0), ...
                      MPB_b*(n_b_1-n_b_0), ...
                      MPB_c*(n_c_1-n_c_0)];
Indirect_Effects_c = [MPC_a*(n_a_1-n_a_0), ...
                      MPC_b*(n_b_1-n_b_0), ...
                      MPC_c*(n_c_1-n_c_0)];
                  
% -------------------------------------------------- %
% (9) Net Worth Component: Bond
% Given bond positions b_0 in state S_0
% What would be the value in S_Star?
% What would be the value in S_1?
% The difference is the bond term
% -------------------------------------------------- %
% Bond and Capital Positions in S_0
sav_a_orig = S_0(jx_sav_a); 
k_a_orig   = S_0(jx_k_a); 
b_a_orig   = S_0(jx_sav_a) - S_0(jx_k_a)*S_0(jx_q); 
sav_b_orig = S_0(jx_sav_b); 
k_b_orig   = S_0(jx_k_b); 
b_b_orig   = S_0(jx_sav_b) - S_0(jx_k_b)*S_0(jx_q); 
sav_c_orig = S_0(jx_sav_c); 
k_c_orig   = S_0(jx_k_c); 
b_c_orig   = S_0(jx_sav_c) - S_0(jx_k_c)*S_0(jx_q); 
% 
% Would-be bond value in S_Star
bv_a_0 = b_a_orig * S_0(jx_nom_i) / S_Star(jx_infl);
bv_b_0 = b_b_orig * S_0(jx_nom_i) / S_Star(jx_infl);
bv_c_0 = b_c_orig * S_0(jx_nom_i) / S_Star(jx_infl);

% Would-be bond value in S_1
bv_a_1 = b_a_orig * S_0(jx_nom_i) / S_1(jx_infl);
bv_b_1 = b_b_orig * S_0(jx_nom_i) / S_1(jx_infl);
bv_c_1 = b_c_orig * S_0(jx_nom_i) / S_1(jx_infl);

% Bond Components
Bond_Component = [bv_a_1 - bv_a_0, bv_b_1 - bv_b_0, bv_c_1 - bv_c_0];

% -------------------------------------------------- %
% (10) Net Worth Component: Capital, Profit
% -------------------------------------------------- %
Profit_Component = [k_a_orig * (S_1(jx_pi)-S_Star(jx_pi)), ...
                    k_b_orig * (S_1(jx_pi)-S_Star(jx_pi)), ...
                    k_c_orig * (S_1(jx_pi)-S_Star(jx_pi))];

% -------------------------------------------------- %
% (11) Capital Price Component
% -------------------------------------------------- %
Capital_Price_Component = (1-prms.delta) * ...
    [k_a_orig * (S_1(jx_q)-S_Star(jx_q)), ...
     k_b_orig * (S_1(jx_q)-S_Star(jx_q)), ...
     k_c_orig * (S_1(jx_q)-S_Star(jx_q))];

% -------------------------------------------------- %
% (12) Transfer Component
% -------------------------------------------------- %
NetWorth_Change    = [n_a_1 - n_a_0, n_b_1 - n_b_0, n_c_1 - n_c_0];
NetWorth_0         = [n_a_0, n_b_0, n_c_0];
NetWorth_1         = [n_a_1, n_b_1, n_c_1];
Transfer_Component = NetWorth_Change - Capital_Price_Component - Profit_Component - Bond_Component;

% -------------------------------------------------- %
% (13) Compute the separate d(wl)/d(veps) term
% -------------------------------------------------- %                
Labor_Income_0 = S_Star(jx_w) * S_Star(jx_l)*[prms.labor_alloc_a*prms.lmbd_a, ...
     prms.labor_alloc_b*prms.lmbd_b, ...
     prms.labor_alloc_c*prms.lmbd_c];
Labor_Income_1 = S_1(jx_w) * S_1(jx_l)*[prms.labor_alloc_a*prms.lmbd_a, ...
     prms.labor_alloc_b*prms.lmbd_b, ...
     prms.labor_alloc_c*prms.lmbd_c];
Labor_Income_Change = S_1(jx_w) * S_1(jx_l) - S_Star(jx_w) * S_Star(jx_l);                
Labor_Income_Component = Labor_Income_Change * ...
    [prms.labor_alloc_a*prms.lmbd_a, ...
     prms.labor_alloc_b*prms.lmbd_b, ...
     prms.labor_alloc_c*prms.lmbd_c];

% -------------------------------------------------- %
% (14) Adjust Direct and Indirect Effect
% Add Labor_Income_Component * MPR to Indirect Effect
% Take them away from direct effect
% -------------------------------------------------- %  
Direct_Effects_k = Direct_Effects_k - Labor_Income_Component .* MP_k;
Direct_Effects_b = Direct_Effects_b - Labor_Income_Component .* MP_b;
Direct_Effects_c = Direct_Effects_c - Labor_Income_Component .* MP_c;

Indirect_Effects_k = Indirect_Effects_k + Labor_Income_Component .* MP_k;
Indirect_Effects_b = Indirect_Effects_b + Labor_Income_Component .* MP_b;
Indirect_Effects_c = Indirect_Effects_c + Labor_Income_Component .* MP_c;

NetWorth_Component = NetWorth_Change + Labor_Income_Component;

% -------------------------------------------------- %
% (15) Compute Residuals
% -------------------------------------------------- %  
Residual_k = Total_Effects_k - Direct_Effects_k - Indirect_Effects_k;
Residual_b = Total_Effects_b - Direct_Effects_b - Indirect_Effects_b;
Residual_c = Total_Effects_c - Direct_Effects_c - Indirect_Effects_c;


Other_Effects_c = Direct_Effects_c - Return_Effects_c;
Other_Effects_b = Direct_Effects_b - Return_Effects_b;
Other_Effects_k = Direct_Effects_k - Return_Effects_k;

%% 6. Further Decompose Return Effect on Capital
% -------------------------------------------------- %
% E(rk) Change
% -------------------------------------------------- %  
Erk_change = (S_1(jx_Erk) - S_Star(jx_Erk));
Erk_component_k(1) = Erk_change * S_Star(jx_dk_dErk1);
Erk_component_k(2) = Erk_change * S_Star(jx_dk_dErk2);
Erk_component_k(3) = Erk_change * S_Star(jx_dk_dErk3);

% -------------------------------------------------- %
% E(rf) Change
% -------------------------------------------------- %  
Erf_change = (S_1(jx_Erf) - S_Star(jx_Erf));
Erf_component_k(1) = Erf_change * S_Star(jx_dk_dErf1);
Erf_component_k(2) = Erf_change * S_Star(jx_dk_dErf2);
Erf_component_k(3) = Erf_change * S_Star(jx_dk_dErf3);


 decomp_mat(:,1,ccc)  = ones(3,1).*irf_factor_1*log(sum(Total_Effects_k + [k_a_0, k_b_0, k_c_0])/sum([k_a_0, k_b_0, k_c_0]));
 decomp_mat(:,2,ccc)  = irf_factor*log(([k_a_0, k_b_0, k_c_0] + Total_Effects_k)./[k_a_0, k_b_0, k_c_0]);
 decomp_mat(:,3,ccc)  = NetWorth_0./(S_0([jx_sav_a, jx_sav_b, jx_sav_c]).*(1-S_0([jx_sh_a, jx_sh_b, jx_sh_c])));
 decomp_mat(:,4,ccc)  = MP_k;
 decomp_mat(:,5,ccc)  = irf_factor*log((NetWorth_0 + Labor_Income_0 + NetWorth_Component)./(NetWorth_0 + Labor_Income_0));
 decomp_mat(:,6,ccc)  = irf_factor*log((Direct_Effects_k+[k_a_0, k_b_0, k_c_0])./[k_a_0, k_b_0, k_c_0]);
 decomp_mat(:,7,ccc)  = decomp_mat(:,2,ccc) - decomp_mat(:,3,ccc).*decomp_mat(:,4,ccc).*decomp_mat(:,5,ccc) - decomp_mat(:,6,ccc);
  
 decomp_mat2(1,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(NetWorth_1(1)/NetWorth_0(1)-sum(NetWorth_1)/sum(NetWorth_0)); 
 decomp_mat2(2,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(Labor_Income_Component(1)/NetWorth_0(1) - sum(Labor_Income_Component)/sum(NetWorth_0));
 decomp_mat2(3,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(Bond_Component(1)/NetWorth_0(1) - sum(Bond_Component)/sum(NetWorth_0));
 decomp_mat2(4,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(Profit_Component(1)/NetWorth_0(1) - sum(Profit_Component)/sum(NetWorth_0));
 decomp_mat2(5,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(Capital_Price_Component(1)/NetWorth_0(1) - sum(Capital_Price_Component)/sum(NetWorth_0));
 decomp_mat2(6,ccc)  = irf_factor_1*NetWorth_0(1)/sum(NetWorth_0)*(Transfer_Component(1)/NetWorth_0(1) - sum(Transfer_Component)/sum(NetWorth_0));
