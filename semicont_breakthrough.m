function t_bt_analytic = semicont_breakthrough(K,Q,n_el,suffix,cin,v_total,t_total,eps,A_mem,nu,module_len,c0,clim,module_mass,rho_mem_SI,dname)
         
   
    % flag for plotting
    plots = false;

    %% Analytic breakthrough time calculation
    % Ref:
    % https://www.annualreviews.org/doi/full/10.1146/annurev-chembioeng-061312-103318
    % OR doi: 10.1146/annurev-chembioeng-061312-103318
    % page 125, paragraph just before Constant pattern section

    % intermediate 1
    na_by_ca = K * Q / (1+K*cin);
    
    % membrane volume
    v_mem = module_mass / rho_mem_SI; %[m3]
    
    % bed volumes
    N_BV = v_mem/v_total;
    
    % Breakthrough time for membrane of defined thickness
    t_bt_analytic = eps * t_total * N_BV * (1 + ((1-eps)/(eps))*((K*Q)/(1+K*cin))); %[s]
    disp(['Analytic breakthrough time: ',num2str(t_bt_analytic/(3600*24)),' days']);

    %% Numeric breakthrough time calculation
    % simulation parameters
    length_span = linspace(0, module_len, n_el);

    del_z = length_span(2) - length_span(1);

    % initial concentration in finite elements
    cinit = c0*ones(n_el,1);

    ode = @(t,c) semicont_ode(t,c,v_total,t_total,eps,A_mem,nu,K,Q,cin,n_el,del_z);

    % integration timespan
    tspan = [0 10*365*24*3600]; % [years]

    % run simulation
    [T C] = ode23(ode,tspan,cinit);

    % change code here to change breakthrough criteria
    cbreak = cin;
    
%     if sum(C(:,end) < cbreak) == 0
    if C(:,end) < cbreak
        disp(' ')
        disp('WARNING: Breakthrough did not happen. Increase simulation horizon.')
        disp([num2str(tspan(end)/3600/24),' days is not long enough.'])
        disp(['Outlet concentration at final time: ',num2str(C(end,end))])
        disp(['Breakthrough concentration criteria: ',num2str(cbreak)])
        disp(' ')
    end
    
    over_limit = C(:,end) > cbreak; 
    %     over_limit = C(end,:) > clim;
    i = find(over_limit, 1);
    disp(['Numeric breakthrough time: ',num2str(T(i)/(3600*24)),' days.'])

    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/semibatch_op_T_',suffix,'.csv'),T);
    csvwrite(strcat(dname,'/semibatch_op_C_',suffix,'.csv'),C);
    csvwrite(strcat(dname,'/semibatch_op_lengths_',suffix,'.csv'),length_span);

    % figure;
    % contour(C)
    if plots
        figure;
        [concs,h] = contour(T/(3600*24),length_span*1e6,C.');
        xlabel('time(days)');
        ylabel('membrane thickness(microns)');
        xlim([0,1.1*T(i)/(3600*24)])
        grid on;

        figure;
        plot(T/(3600*24),C(:,end)*1e3);
        xlabel('time(days)');
        ylabel('fluid phase concentration(mmol/cm3)');
    end
end