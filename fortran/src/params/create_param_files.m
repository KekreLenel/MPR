% -------------------------------------------------------------------------
% create_param_files.m: creates parameter files for comparative calibrations  
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

% Benchmark parameters
% ======================================================
% ======================================================

% Population Shares 
prm.lmbd_a          = 0.04; % pop share a
prm.lmbd_b          = 0.36; % pop share b

% Preference Parameters
prm.bbeta_a         = 0.98;
prm.bbeta_b         = 0.98;
prm.bbeta_c         = 0.98;
prm.gma_a           = 10.0; 
prm.gma_b           = 25.5;
prm.gma_c           = (prm.lmbd_a/(prm.lmbd_a + prm.lmbd_b)/prm.gma_a + prm.lmbd_b/(prm.lmbd_a + prm.lmbd_b)/prm.gma_b)^(-1);
prm.ies_a           = 0.8;
prm.ies_b           = 0.8;
prm.ies_c           = 0.8;
prm.theta           = 1.0;  % frisch     

% Technology Parameters
prm.delta           = 0.025;
prm.aalpha          = 0.33;
prm.sig_z           = 0.0055;   
prm.chiX            = 3.5;      % capital adjustment cost

% Monetary Policy
prm.phi             = 1.5;      % Taylor rule coefficient
prm.tayl_ic         = 0.003;    % Taylor rule intercept: approx. 0% inflation 
prm.sig_m           = 0.0025/4; % Monetary policy shock
prm.rho_m           = 0.0;

% Disaster Shock
disast_p            = 0.005;      % average disaster probability 
disast_std          = 0.0025;     % unconditional standard deviation 
[X, FVAL, FLAG]     = fsolve(@(x) 1000*[ disast_p - exp(x(1) + x(2)/2) ; disast_std^2 - (exp(x(2))-1)*exp(2*x(1) + x(2)) ], [-1; 1]);
prm.disast_p        = X(1);       % average log(p)
prm.varphi_p        = 0.15;       % disaster depth
prm.rho_p           = 0.8;        % disaster autocorrelation
prm.disast_std      = sqrt(X(2)); % unconditional std. log(p)

% Rotemberg wage adjustment
prm.vareps_w        = 10;
prm.tau_w           = 1/(1-prm.vareps_w);
prm.chiW            = 150;

% Newborn wealth allocation 
prm.s_bar_a         =  0.0;
prm.s_bar_c         = -0.25; 

% Target wealth shares (a and c wealth are state variables)
prm.s_trgt_a        = 0.18;
prm.s_trgt_c        = 0.23;

% Mortality rate
prm.xi              = 0.01;

% Average labor supply 
prm.l_target        = 1.0;

% Labor share: phi_i
prm.labor_alloc_a   = 0.03/prm.lmbd_a; 
prm.labor_alloc_b   = 0.14/prm.lmbd_b; 
prm.labor_alloc_c   = 0.83/(1.0 - prm.lmbd_a - prm.lmbd_b); 

% lower bound capital holdings
prm.kbar            =  10;

% Government debt 
prm.gov_debt        = 0.10;

% Grid scaling parameters: width and center relative to steady state
prm.k_grid_adj      = 1.25;
prm.w_grid_adj      = 1.06;
prm.s_a_dev         = 0.1;
prm.s_c_dev         = 0.1;
prm.k_dev_param     = 0.1;
prm.w_dev_param     = 0.05; 

% IRF shock sizes 
prm.IRF_g           =  2.0;
prm.IRF_m           = -10.0;
prm.IRF_dis         =  2.0;

% for efficiency kbar constrained only activated for group c, not binding for a anddkj b 
prm.constrained_a   = 0.0;
prm.constrained_b   = 0.0;
prm.constrained_c   = 1.0;

% idiosyncratic risk
prm.use_idio_risk   = 0.00;
prm.idio_risk_a     = 0.00;
prm.idio_risk_b     = 0.00;
prm.idio_risk_c     = 0.00;

% ======================================================

% number of parameters 
n_params = size(fieldnames(prm),1);

% how many alternative calibrations 
n_comp = 9;

% store number of calibrations in file
fid = fopen(['..', filesep, 'output', filesep, 'tmp', filesep, 'n_comp.txt'], 'w');
fprintf(fid, '%i', n_comp);
fclose(fid);

% how many variables varied in comp. statics
n_comp_prms = 33;

% Which parameters to vary
prm_list = {'rho_m', 'chiX', 'chiW', 'use_idio_risk', 'idio_risk_a','idio_risk_b', ...
            'idio_risk_c', 'bbeta_a', 'bbeta_b', 'bbeta_c', 'gma_a', 'gma_b', 'gma_c', ...
            'constrained_c', 'xi', 's_bar_a', 's_bar_c', 'lmbd_a', 'lmbd_b', ...
            's_trgt_a', 's_trgt_c', 'labor_alloc_a', 'labor_alloc_b', ...
            'labor_alloc_c', 'k_grid_adj', 'w_grid_adj', 's_a_dev', 's_c_dev',  ... 
            'tayl_ic', 'disast_p', 'varphi_p', 'rho_p', 'disast_std'};

comp_mat(1:n_comp, 1)  = prm.rho_m;
comp_mat(1:n_comp, 2)  = prm.chiX;
comp_mat(1:n_comp, 3)  = prm.chiW;
comp_mat(1:n_comp, 4)  = prm.use_idio_risk;
comp_mat(1:n_comp, 5)  = prm.idio_risk_a;
comp_mat(1:n_comp, 6)  = prm.idio_risk_b;
comp_mat(1:n_comp, 7)  = prm.idio_risk_c;
comp_mat(1:n_comp, 8)  = prm.bbeta_a;
comp_mat(1:n_comp, 9)  = prm.bbeta_b;
comp_mat(1:n_comp, 10) = prm.bbeta_c;
comp_mat(1:n_comp, 11) = prm.gma_a;
comp_mat(1:n_comp, 12) = prm.gma_b;
comp_mat(1:n_comp, 13) = prm.gma_c;
comp_mat(1:n_comp, 14) = prm.constrained_c;
comp_mat(1:n_comp, 15) = prm.xi;
comp_mat(1:n_comp, 16) = prm.s_bar_a;
comp_mat(1:n_comp, 17) = prm.s_bar_c;
comp_mat(1:n_comp, 18) = prm.lmbd_a;
comp_mat(1:n_comp, 19) = prm.lmbd_b;
comp_mat(1:n_comp, 20) = prm.s_trgt_a;
comp_mat(1:n_comp, 21) = prm.s_trgt_c;
comp_mat(1:n_comp, 22) = prm.labor_alloc_a;
comp_mat(1:n_comp, 23) = prm.labor_alloc_b;
comp_mat(1:n_comp, 24) = prm.labor_alloc_c;
comp_mat(1:n_comp, 25) = prm.k_grid_adj;
comp_mat(1:n_comp, 26) = prm.w_grid_adj;
comp_mat(1:n_comp, 27) = prm.s_a_dev;
comp_mat(1:n_comp, 28) = prm.s_c_dev;
comp_mat(1:n_comp, 29) = prm.tayl_ic;
comp_mat(1:n_comp, 30) = prm.disast_p;     
comp_mat(1:n_comp, 31) = prm.varphi_p; 
comp_mat(1:n_comp, 32) = prm.rho_p; 
comp_mat(1:n_comp, 33) = prm.disast_std;

% (2) RANK
comp_mat(2, 11:13) = (prm.s_trgt_a/prm.gma_a  + (1 - prm.s_trgt_a - prm.s_trgt_c)/prm.gma_b+ prm.s_trgt_c/prm.gma_c)^(-1); % risk aversion 
comp_mat(2, 14)    = 0.0;        % no constraints 
comp_mat(2, 16:21) = 1/3;        % lmbd, s
comp_mat(2, 22:24) = 1.0;        % labor_alloc 
comp_mat(2, 25)    = 1.05;       % k_grid_adj 
comp_mat(2, 26)    = 1.02;       % w_grid_adj  
comp_mat(2, 29)    = 0.005;      % tayl_ic 

% (2) MP Shock Persistence
comp_mat(3, 1)  = 0.75;

% (3) No adjustment cost
comp_mat(4, 2) = 0.0;

% (4) No wage rigidity
comp_mat(5, 3) = 0;

% (6) Idiosyncratic risk 
comp_mat(6,4)     = 1.0;    % use idiosyncratic risk 
comp_mat(6,5)     = 0.0;    % idio on a
comp_mat(6,6)     = 0.036;  % idio on b        
comp_mat(6,7)     = sqrt(prm.lmbd_a/(prm.lmbd_a + prm.lmbd_b)*comp_mat(6,5) + prm.lmbd_b/(prm.lmbd_a + prm.lmbd_b)*(comp_mat(6,6)^2));  % idiosyncratic risk on c 
comp_mat(6,8:10)  = 0.98;  % bbeta 
comp_mat(6,11:13) = 11.0;   % risk aversion 
comp_mat(6,15)    = 0.01;   % xi
comp_mat(6,16)    = -0.02;  % s_bar_a
comp_mat(6,17)    = -0.15;  % s_bar_c
comp_mat(6,25)    = 1.20;   % k_grid_adj 
comp_mat(6,26)    = 1.05;   % w_grid_adj  
comp_mat(6,29)    = 0.0035; 
disast_p          = 0.005;
disast_std        = 0.009;
[X, FVAL, FLAG]   = fsolve(@(x) 1000*[ disast_p - exp(x(1) + x(2)/2) ; disast_std^2 - (exp(x(2))-1)*exp(2*x(1) + x(2)) ], [-1; 1]);
comp_mat(6,30)    = X(1);
comp_mat(6,33)    = sqrt(X(2));

% (7) Idiosyncratic risk (rep agent)
comp_mat(7, 4)     = 1.0;     % use idio risk 
comp_mat(7, 5:7)   = (prm.s_trgt_a*comp_mat(6,5) +  (1-prm.s_trgt_a - prm.s_trgt_c)*comp_mat(6,6)+ prm.s_trgt_c*comp_mat(6,7));          
comp_mat(7, 8:10)  = comp_mat(6,8);  % beta
comp_mat(7, 11:13) = comp_mat(6,11); % gma
comp_mat(7, 15)    = comp_mat(6,15); % xi
comp_mat(7, 14)    = 0.0;     % constraint  
comp_mat(7, 16:21) = 1/3;     % lmbd, s
comp_mat(7, 22:24) = 1.0;     % labor_alloc 
comp_mat(7, 25)    = 1.05;    % k_grid_adj 
comp_mat(7, 26)    = 1.02;    % w_grid_adj  
comp_mat(7, 29)    = 0.006;  % tayl_ic 
comp_mat(7, 30)     = X(1);
comp_mat(7, 33)     = sqrt(X(2));

% hedge fund calibration 
comp_mat(8, 8:10)  = 0.984; % bbeta 
comp_mat(8, 11)    = 2.5;  % risk aversion 
comp_mat(8, 12)    = 21.0; % risk aversion 
comp_mat(8, 13)    = (0.004/0.4/2.5 + 0.396/0.4/21)^(-1);  
comp_mat(8, 15)    = 0.04;  %  xi
comp_mat(8, 16)    = 0.0;  % s_bar_a
comp_mat(8, 17)    = 0.009/0.04;  % s bar_c
comp_mat(8, 18)    = 0.004; % lmbd_a
comp_mat(8, 19)    = 0.396; % lmbd_b
comp_mat(8, 20)    = 0.02;  % s_trg_a 
comp_mat(8, 21)    = 0.23;  % s trg_c
comp_mat(8, 22)    = 0.00/0.004;  % labor alloc 
comp_mat(8, 23)    = 0.17/0.396;  % 
comp_mat(8, 24)    = 0.83/0.6;    %
comp_mat(8, 25)    = 1.08;   % k_grid_adj 
comp_mat(8, 26)    = 1.02;   % w_grid_adj  
comp_mat(8, 27)    = 0.019; % s_a_dev
% comp_mat(8, 28)    = 0.08;  % s_c_dev
comp_mat(8, 29)    = 0.003; % tayl_ic

% hedge fund calibration (rep)
comp_mat(9, 8:10)  = comp_mat(8,8);
comp_mat(9, 11:13) = (comp_mat(8, 20)/comp_mat(8,11) +  (1 - comp_mat(8, 20) - comp_mat(8, 21))/comp_mat(8,12)+ comp_mat(8, 21)/comp_mat(8,13))^(-1); %19.0;  % risk aversion 
comp_mat(9, 14)    = 0.0;   % constraint type 
comp_mat(9, 16:21) = 1/3;        % lmbd, s
comp_mat(9, 22:24) = 1.0;        % labor_alloc 
comp_mat(9, 25)    = 1.05;  % k_grid_adj 
comp_mat(9, 26)    = 1.02;  % w_grid_adj  
comp_mat(9, 29)    = 0.006; % tayl_ic 
comp_mat(9, 15)    = comp_mat(8,15);


% Now generate parameter files
for nnn = 1:n_comp

    name_list = string(fieldnames(prm));

    for ppp = 1:n_comp_prms
        prm = setfield(prm, prm_list{ppp}, comp_mat(nnn,ppp) );
    end

    prm_vec = cell2mat(struct2cell(prm));

    file_str = ['..', filesep, 'src', filesep, 'params', filesep, sprintf('param_file_%i.csv',nnn)];
    fid = fopen(file_str, 'w');
    fprintf(fid, [ repmat('%s, ', [1, n_params-1]  ), ' %s\n'], name_list(1:end-1), name_list(end));
    fprintf(fid, [ repmat('%10.6f, ', [1, n_params-1]  ), ' %10.6f\n'], prm_vec(1:end-1), prm_vec(end));
    fclose(fid);

end

