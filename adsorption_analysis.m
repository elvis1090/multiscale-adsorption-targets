function adsorption_analysis(v_total,t_total,t_regen,module_mass,dname)
% This function runs all the analysis for adsorpive membrane process / paper
    %% Input data, units indicated within []

    N = 1.6e14; % [1/m2] number of pores per m2 of membrane
    dp = 5.1e-9; % [m] pore diameter
    % N and dp data taken from supplementary information of DOI: 10.1039/C6TA09287J 

    % viscosity of water
    mu = 8.9e-4; % [Pa-s] density of water at room temperature

    % density of polymer solution
    rho_polym = 1; % [g/cm3] Approximate average density of a polymer solution
    % since polyisoprene: 0.92 g/mL  https://polymerdatabase.com/polymers/polyisoprene.html
    % and polystyrene: 1.05 g/mL https://polymerdatabase.com/polymers/polystyrene.html

    N_A = 6.022e20; %[1/mmol] Avogadro's number for 1 mmol of substance

    mw_Pb = 207.2e-3; %[g/mmol] Molecular mass of 1 mmol of lead i.e. target solute / contaminant

    rho_Pb = 11.34; %[g/cm3] Density of lead i.e. target solute / contaminant

    % save water_consumption 
    csvwrite(strcat(dname,'/v_total.csv'),v_total);

    Am_old = 50.71; % [m2] membrane area, area of 1 kg of membrane which was the mass used in the batch analysis.
    A_mem = ((pi/4)*(2.54e-2)^2)/(17.6e-6); % [m2] membrane area for 1 kg of membrane
    % Area (A_mem) calculated abased on information that 1" dia. of membrane weighs
    % 17.6 mg found in SI of doi: 10.1021/acs.langmuir.5b01605

    l_mem = 50e-6; % [m] membrane thickness
    % l_mem taken from value reported in the SI of DOI: 10.1021/acsami.7b04603

    act_thk = 40e-9; % [m] active layer thickness 
    % act_thk taken from value reported in the SI of DOI: 10.1039/c6ta09287j

    del_P_psi = 65; % [psi] pressure drop for a conventional household RO membrane
    % ref: spec sheet for pentair tlc-24/36/50/75/100 ro membrane filter in DowlingLab/simple-batch/current_tech/

    cin_ppb = 100; %[ppb] Lead concentration in Flint Michigan water 
    % ref: 1) https://www.mlive.com/news/flint/2016/03/some_flint_water_test_sites_st.html
    % 2) doi: 10.1021/acsami.7b04603

    c0 = 1e-9; %[mmol/cm3] initial bed concentration. 0 for a fully regenerated bed

    clim_ppb = 15; %[ppb] EPA limits for lead concentration in drinking water
    % Ref: https://www.epa.gov/ground-water-and-drinking-water/national-primary-drinking-water-regulations#seven
    % accessed 28-Aug-2019

    sp_thr = 1e-3; %[m] The maximum (threshold) value of membrane thickness for it to be formed into 
    % a spiral wound module

    % arrays of K and Q values to generate breakthrough curves
    K_paper_values = [0.87 6400*1e-3]; % [l/mmol] binding affinity 
    % Ref: DOI: 10.1021/acsami.7b04603, 10.1021/acscentsci.8b00690 

    Q_paper_values = [1.38 1.2]; % [mmol/g] saturation capacity 
    % Ref: DOI: 10.1021/acsami.7b04603, 10.1021/acscentsci.8b00690

    % array of filename suffixes to save breakthrough data, corresponding to
    % respective K and Q values
    suffixes = ["pash" "terp"];

    % masses for material property target calculations 
    mpt_masses = [1 10 25 50 100 1000]; % [kg] 

    % number of elements for curve calculations
    n_el_mpt = 100000;
    
    % number of elements for breakthrough integration
    n_el_bt = 1000;

    % linear version of mpt_masses to calculate curve bounding saturation
    % capacity (Q)
    cont_mpt_masses = logspace(log10(mpt_masses(1)),log10(mpt_masses(end)*40),n_el_mpt);
    % cont_mpt_masses = linspace(mpt_masses(1),mpt_masses(end)*10,n_el_mpt);

    % pressure set points for household water supply for capacity upper bound
    % calculations
    % ref:
    % https://www.plumbingsupply.com/residential-water-pressure-explained.html,
    % accessed 22-jan-2020
    % [upper limit, approx middle, average]
    del_P_vec = [45 65 80]; % psi

    % save to file for plotting
    writematrix(del_P_vec,strcat(dname,'/del_P_vec.csv'));

    % masses for pore-size influence calculations
    ps_masses = linspace(1,1000,1000);

    % range of K for mpt_curve calculations
    % K_mpt = [linspace(1e-1, 1e3, n_el_mpt/2) linspace(1e3, 1e8, n_el_mpt/2)]; %[cm3/mmol = l/mol = 1/M]    
    K_mpt = logspace(-1, 8, n_el_mpt); %[cm3/mmol = l/mol = 1/M]    


    % parameters for sensitivity analysis over dimensionless material
    % property target curves
    % rec_span - values of recovery ratio for the sensitivity analysis, say
    % [r1, r2, r3]
    % eps_span - values of porosity for the sensitivity analysis,
    % corresponding to rec_span, say [e1, e2, e3]
    % N_BV_span_all - values of bed volume (N_BV) to print on a plot for a
    % given a combination of recovery ratio and porosity. Say N_BV_span_all
    % = [nbv_1a, nbv_1b, nbv_1c; nbv_2a, nbv_2b, nbv_2c; nbv_3a, nbv_3b, nbv_3c]
    % sensitivity analysis will provide data for 3 plots. 
    % plot 1: 3 lines for nbv_1a, nbv_1b, nbv_1c at fixed r1 and e1
    % plot 2: 3 lines for nbv_2a, nbv_2b, nbv_2c at fixed r2 and e2
    % plot 3: 3 lines for nbv_3a, nbv_3b, nbv_3c at fixed r3 and e3
    % plots are generated using associated python script
    rec_span = [1.1 10 100 14]; %[-]
    eps_span = [0.3, 0.3, 0.3, 0.3]; %[-]
%     rec_span = [1.1 2.5 100 14]; %[-]
%     eps_span = [0.3, 0.3, 0.3, 0.3]; %[-]
    N_BV_span_all = [10, 500, 50000; 10, 500, 50000; 5, 15, 25; 1000, 5000, 10000];
    writematrix(rec_span,strcat(dname,'/rec_span.csv'));
    writematrix(N_BV_span_all,strcat(dname,'/N_BV_span_all.csv'));
    writematrix(eps_span,strcat(dname,'/eps_span.csv'));
    
    % Vector to store feasible limit of Q calculations
    Q_feaslim = zeros(length(mpt_masses),1);
    Q_feaslim3 = zeros(length(mpt_masses),1);

    %% model parameter calculations
    disp('');
    disp('model parameter calculations');

    % void fraction [-] = pore area / membrane area
    % eps = N * pi * dp^2 / (4*1);
    % Note: hard coding void fraction / porosity based on a value Bill
    % mentioned during discussions on 03-Oct-2019. See Research meeting notes
    % on 04-Oct-2019 for more details.
    eps = 0.3;

    % phase ratio [-]
    nu = (1-eps)/eps;

    % intersitial velocity [m/s] = flow_rate / (bed void fraction * membrane area)
    % V = water_flow / (eps*A_mem); %[m/s]
    % will move to functions where it is being used.

    % K [cm3/mmol]
    K = K_paper_values .* 1e3;  %[cm3/mmol]

    % Area of 1 inch diameter of membrane 
    A_mem_1 = (pi/4)*(2.54e-2)^2; % [m2]

    % Volume of 1 inch diameter of membrane
    vol_mem_1 = A_mem_1 * l_mem; %[m3]

    % density of the membrane
    rho_mem = (17.6e-3)/(vol_mem_1*(1e2)^3); %[g/cm3]

    % density calculations based on information that 1" dia. of membrane weighs
    % 17.6 mg found in SI of doi: 10.1021/acs.langmuir.5b01605 and l_mem taken 
    % from value reported in the SI of DOI: 10.1021/acsami.7b04603

    % density of the membrane in SI units
    rho_mem_SI = rho_mem * 1e-3 / (1e-2)^3; %[kg/m3]
    % print rho_mem_SI to text file for use in plotting notebook
    csvwrite(strcat(dname,'/rho_mem_SI.csv'),rho_mem_SI);

    % volume of a mem_module 
    module_vol = module_mass / rho_mem_SI;

    % module length
    module_len = module_vol / A_mem;

    % Q [mmol/cm3]
    Q = Q_paper_values .* rho_mem; %[mmol/cm3]

    % cin [mmol/cm3]
    cin = cin_ppb * (1/207.2) / (1e3); %[mmol/cm3]
    clim = clim_ppb * (1/207.2) / (1e3); %[mmol/cm3]

    % dimensionless K [-]
    K_nondim = logspace(-4,4,n_el_mpt);

    %% Sensitivity analysis over non-dimensional MPT curves
    % This is a heavy calculation and will be performed only if case4 is
    % being run
    
    if strcmp(dname,'case4')          
        %% sensitivity analysis over MPT curves
        
        disp('Sensitivity analysis over non-dimensional MPT curves');
        
        %save K_nondim to a .csv file to plot in Python
        csvwrite(strcat(dname,'/sbmpt_Knondim.csv'),K_nondim.');
        
        for i = 1:length(rec_span)
            N_BV_span = N_BV_span_all(i,:);
            
            for j = 1:length(N_BV_span)
                nondim_sens(eps_span(i),N_BV_span(j),rec_span(i),K_nondim,1,dname)                
            end
        end
        
        %% Fit Langmuir isotherm for arsenic adsorbent HAX1 
        % Nonlinear regression for isotherm
        disp('Nonlinear regression for HAX1 Isotherm');
        
        % test change read isotherm data and extract into meaningful variables
        % ref: Fig 1C of paper with doi:
        % http://dx.doi.org/10.1016/j.scitotenv.2017.05.126
        dat_hax1 = readmatrix('arsenic_data_fisotherm.csv');      
        c_iso_hax1 = dat_hax1(:,1); % mM
        q_iso_hax1 = dat_hax1(:,2); %[mmol/g_membrane]
        
        q_pass = q_iso_hax1;
        c_pass = c_iso_hax1;
        NL_reg = @(z) isotherm_regression(z,q_pass,c_pass);

        %Initial guesses for Q and K are 0
        guess = [0 0];

        %Fit nonlinear regression parameters using fminsearch(), i.e. solve an 
        %unconstrained nonlinear optimization problem.
        z = fminsearch(NL_reg,guess);

        %Save and display results
        Q_NL = z(2);
        K_NL = z(1);
        disp('HAX Nonlinear Model:');
        disp(['Q=',num2str(Q_NL),'(mmol/g_membrane)']);
        disp(['K=',num2str(K_NL),'(l/mmol)']);
        disp(' ');

        csvwrite(strcat(dname,'/Q_fit_hax1.csv'),Q_NL);
        csvwrite(strcat(dname,'/K_fit_hax1.csv'),K_NL);
        csvwrite(strcat(dname,'/c_iso_hax1.csv'),c_iso_hax1);
        csvwrite(strcat(dname,'/q_iso_hax1.csv'),q_iso_hax1);
        
        %% Calculate dimensionless properties for arsenic case study
        calc_dmlss_props(dname);
        
    elseif strcmp(dname,'case5')
        %% Analyze over-/under-predictions to due incorrect scaleup
        analyze_overpredictions(c0,v_total,t_total,eps,A_mem,nu,cin,rho_mem_SI,dname)
       
    elseif strcmp(dname,'case6')
        %% Plots for the Supporting Information
        
        %% calculate the breakthrough time for -PASH and -Terp membrane materials
        disp('');
        disp('calculate the breakthrough time for -PASH and -Terp membrane materials');

        % define an anonymous function to pass constant data to the
        % semicont_breakthrough function
        breakthrough = @(K,Q,suff) semicont_breakthrough(K,Q,n_el_bt,suff,cin,v_total,t_total,eps,A_mem,nu,module_len,c0,clim,module_mass,rho_mem_SI,dname);

        for i = 1:length(suffixes)
            disp(suffixes(i));
            notused = breakthrough(K(i), Q(i), suffixes(i));
        end
        
        %% Batch material property target curves
        disp('')
        disp('Batch material property targets');
        bat_mpt = @(mpt_mass) batch_mpt(cin,clim,K_mpt,v_total,mpt_mass,rho_mem,rho_mem_SI,eps,dname);

        for j = 1:length(mpt_masses)
            bat_mpt(mpt_masses(j));
        end
        
        % save to file
        csvwrite(strcat(dname,'/batmpt_Kspan.csv'),K_mpt);
        
        %% Schematics for derivation of mathematical regions of competitive 
        %% batch and semi-continuous behavior
        
        % parameters for sensitivity analysis over dimensionless material
        % property target curves (see earlier for details)
        rec_span = [1.1 2.5 100 14]; %[-]
        eps_span = [0.3, 0.3, 0.3, 0.3]; %[-]
        N_BV_span_all = [10, 500, 50000; 10, 500, 50000; 5, 15, 25; 1000, 5000, 10000];
        writematrix(rec_span,strcat(dname,'/rec_span.csv'));
        writematrix(N_BV_span_all,strcat(dname,'/N_BV_span_all.csv'));
        writematrix(eps_span,strcat(dname,'/eps_span.csv'));
        
        disp('');
        disp('Schematics for competitive behavior derivation');
        
        %save K_nondim to a .csv file to plot in Python
        csvwrite(strcat(dname,'/sbmpt_Knondim.csv'),K_nondim.');
        
        for i = 1:length(rec_span)
            N_BV_span = N_BV_span_all(i,:);
            
            for j = 1:length(N_BV_span)
                nondim_sens(eps_span(i),N_BV_span(j),rec_span(i),K_nondim,1,dname)                
            end
        end
    elseif (strcmp(dname,'case7'))
        %% Supporting information part 2: Nonlinear regression for Langmuir isotherm
        % Nonlinear regression for isotherm
        %terp
        disp('Nonlinear regression for Terp Isotherm');
        c_iso_terp = [0.92 0.4725 0.24208 0.03225 0.564].'; % [mM]
        q_iso_terp = [1.05064 0.88103 0.79385 0.2 0.9895].'; %[mmol/g_membrane]
        
        q_pass = q_iso_terp;
        c_pass = c_iso_terp;
        NL_reg = @(z) isotherm_regression(z,q_pass,c_pass);

        %Initial guesses for Q and K are 0
        guess = [0 0];

        %Fit nonlinear regression parameters using fminsearch(), i.e. solve an 
        %unconstrained nonlinear optimization problem.
        z = fminsearch(NL_reg,guess);

        %Save and display results
        Q_NL = z(2);
        K_NL = z(1);
        disp('TERP Nonlinear Model:');
        disp(['Q=',num2str(Q_NL),'(mmol/g_membrane)']);
        disp(['K=',num2str(K_NL),'(l/mmol)']);
        disp(' ');

        csvwrite(strcat(dname,'/Q_fit_terp.csv'),Q_NL);
        csvwrite(strcat(dname,'/K_fit_terp.csv'),K_NL);
        csvwrite(strcat(dname,'/c_iso_terp.csv'),c_iso_terp);
        csvwrite(strcat(dname,'/q_iso_terp.csv'),q_iso_terp);

        %pash
        disp('Nonlinear regression for Pash Isotherm');
        c_iso_pash = [0.4 1.8 3.9 10 20 40].'; % [mM]
        q_iso_pash = [0.44 0.85 0.95 1.31 1.25 1.34].'; %[mmol/g_membrane]
        
        q_pass = q_iso_pash;
        c_pass = c_iso_pash;
        NL_reg = @(z) isotherm_regression(z,q_pass,c_pass);

        %Initial guesses for Q and K are 0
        guess = [0 0];

        %Fit nonlinear regression parameters using fminsearch(), i.e. solve an 
        %unconstrained nonlinear optimization problem.
        z = fminsearch(NL_reg,guess);

        %Save and display results
        Q_NL = z(2);
        K_NL = z(1);
        disp('PASH Nonlinear Model:');
        disp(['Q=',num2str(Q_NL),'(mmol/g_membrane)']);
        disp(['K=',num2str(K_NL),'(l/mmol)']);
        disp(' ');

        csvwrite(strcat(dname,'/Q_fit_pash.csv'),Q_NL);
        csvwrite(strcat(dname,'/K_fit_pash.csv'),K_NL);
        csvwrite(strcat(dname,'/c_iso_pash.csv'),c_iso_pash);
        csvwrite(strcat(dname,'/q_iso_pash.csv'),q_iso_pash);
        
    elseif (strcmp(dname,'case8'))
        %% Supporting information part 3 - several MPT curves
        
        disp('')
        disp('Generating reference MPT curves for SI');
        
        % parameters for sensitivity analysis over dimensionless material
        % property target curves (see earlier for details)
        rec_span = [1.5, 5, 10, 50, 100, 1000,...
                    1.5, 5, 10, 50, 100, 1000,...
                    1.5, 5, 10, 50, 100, 1000,...
                    1.5, 5, 10, 50, 100, 1000]; %[-]
        eps_span = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1,...
                    0.2, 0.2, 0.2, 0.2, 0.2, 0.2,...
                    0.3, 0.3, 0.3, 0.3, 0.3, 0.3,...
                    0.4, 0.4, 0.4, 0.4, 0.4, 0.4]; %[-]
        N_BV_span_all = [50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;...
                         50, 500, 1000, 5000, 10000, 50000;];
                         
        writematrix(rec_span,strcat(dname,'/rec_span.csv'));
        writematrix(N_BV_span_all,strcat(dname,'/N_BV_span_all.csv'));
        writematrix(eps_span,strcat(dname,'/eps_span.csv'));
        
        %save K_nondim to a .csv file to plot in Python
        csvwrite(strcat(dname,'/sbmpt_Knondim.csv'),K_nondim.');
        
        for i = 1:length(rec_span)
            N_BV_span = N_BV_span_all(i,:);
            disp(['rec=',num2str(rec_span(i)),' eps=',num2str(eps_span(i))])
            
            for j = 1:length(N_BV_span)
                nondim_sens(eps_span(i),N_BV_span(j),rec_span(i),K_nondim,1,dname)                
            end
        end
        
        
    else
       
        %% calculate Semi-continuous Materials Property (K & Q) Target Curves for a fixed mass of membrane material
        disp('');
        disp('Materials Property (K & Q) Target Curves for a fixed mass of membrane material');

        % semi-continuous process
        sb_mpt = @(mpt_mass) semicont_mpt(rho_mem_SI,mpt_mass,cin,rho_mem,K_mpt,v_total,t_total,t_regen,eps,dname);

        for j = 1:length(mpt_masses)
            sb_mpt(mpt_masses(j));
        end

        % save cin, mpt_masses, and K_mpt to a .csv file (to plot in python)
        csvwrite(strcat(dname,'/sbmpt_cin.csv'),cin);
        csvwrite(strcat(dname,'/sbmpt_mspan.csv'),mpt_masses);      
    
        %% Pore size calculations
        del_P = del_P_psi * 6894.76; % [Pa]

        % minimum pore size required to achieve target separation
        dp_min_SI = semicont_dp(v_total,t_total,mu,module_mass,A_mem,del_P,N,rho_mem_SI); % [m]
        dp_min = dp_min_SI * 1e2;

        %% Capacity calculation
        % maximum possible membrane capacity for target separation

        % molar volume of Pb [cm3/mmol]
        v_bar = mw_Pb/rho_Pb;

        Q_max = semicont_maxcapacity(eps,rho_polym,dp_min,N_A,v_bar);

        %% Sensitivity analysis on maximum membrane capacity
        disp('');
        disp('Sensitivity analysis on maximum membrane capacity');

        % iteration limit
        Qmax_iter = length(ps_masses);

        % vectors to save data
        Qmax_sens = zeros(Qmax_iter);
        dp_min_sens = zeros(Qmax_iter);

        for i = 1:Qmax_iter
            dp_min_sens(i) = semicont_dp(v_total,t_total,mu,ps_masses(i),A_mem,del_P,N,rho_mem_SI); %[m]
            Qmax_sens(i) = semicont_maxcapacity(eps,rho_polym,dp_min_sens(i)*1e2,N_A,v_bar); %[mmol/g]
        end

        %% Sensitvity analysis on pressure drop across membrane
        disp('');
        disp('Sensitvity analysis on pressure drop across membrane');

        % Pressure drop range
        del_P_span = linspace(2.8e5,5.6e5,1000); %[Pa]

        % Calculate pore size requirements for varying pressure drops in the
        % semi-continuous process

        % empty matrix to store results
        % TESTING ONLY
        dptargets_save = zeros(length(mpt_masses),length(del_P_span));

        % declare an anonymous function to calculate pore diameter (dp) targets
        sb_dpt = @(m_mem) semicont_dp_sens(v_total,t_total,mu,m_mem,A_mem,del_P_span,N,rho_mem_SI,dname);

        % calculate pore diameter targets
        for i = 1:length(mpt_masses)
            sb_dpt(mpt_masses(i));
        end

        % save del_P_span to a .csv file (to plot in python)
        csvwrite(strcat(dname,'/sbdpt_delP.csv'),del_P_span);

        %% Sensitvity analysis on pressure drop across membrane for parallel membrane modules
        disp('');
        disp('Sensitvity analysis on pressure drop across membrane for parallel membrane modules');

        % declare an anonymous function to calculate pore diameter (dp) targets in
        % parallel modules
        sb_dpt_parallel = @(m_mem) semicont_dp_parallel_sens(v_total,t_total,mu,m_mem,A_mem,del_P_span,N,rho_mem_SI,sp_thr,dname);

        % calculate pore diameter targets for parallel modules
        for i = 1:length(mpt_masses)
            sb_dpt_parallel(mpt_masses(i));
        end

        %% Feasible Q region calculation - 65 psig
        disp('Feasible Q region calculation');

        % loop over masses being considered
        for i = 1:length(mpt_masses)
            dp_min_SI2 = semicont_dp_parallel(v_total,t_total,mu,mpt_masses(i),A_mem,del_P,N,rho_mem_SI,sp_thr); % [m]
            dp_min2 = dp_min_SI2 * 1e2; % [cm]
            Q_feaslim(i) = semicont_maxcapacity(eps,rho_polym,dp_min2,N_A,v_bar); % [mmol/g]
        end

        % save Q feasible limits list to file
        csvwrite(strcat(dname,'/Q_feaslim.csv'),Q_feaslim);

        % calculate the curve showing the upper bound on Q based on pressure drop
        % (and pore size) considerations

        % declare an anonymous function to pass data
        qub_curve = @(press_drp) semicont_capacity_ub_curve(v_total,t_total,t_regen,mu,cont_mpt_masses,A_mem,press_drp,N,rho_mem_SI,sp_thr,eps,rho_polym,N_A,v_bar,cin,rho_mem,dname);

        % loop over pressures considered
        for i=1:length(del_P_vec)
            [K_ub, Q_ub] = qub_curve(del_P_vec(i)); 
        end

        % 80 psig
        del_P3 = 80 * 6894.76; % [Pa]

        % loop over masses being considered
        for i = 1:length(mpt_masses)
            dp_min_SI3 = semicont_dp_parallel(v_total,t_total,mu,mpt_masses(i),A_mem,del_P3,N,rho_mem_SI,sp_thr); % [m]
            dp_min3 = dp_min_SI3 * 1e2; % [cm]
            Q_feaslim3(i) = semicont_maxcapacity(eps,rho_polym,dp_min3,N_A,v_bar); % [mmol/g]
        end

        % save Q feasible limits list to file
        csvwrite(strcat(dname,'/Q_feaslim3.csv'),Q_feaslim3);
        
        %% Calculate masses of adsorbents required for the process
        disp('Calculate existing adsorbent mass requirement')
        calc_ads_masses(v_total, cin, eps, rho_mem, clim, dname);
        
    end % analysis of cases 1 to 3
end % function adsorption_analysis