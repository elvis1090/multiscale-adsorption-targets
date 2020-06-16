function batch_mpt(cin,cout,K_span,v_total,mpt_mass,rho_mem,rho_mem_SI,eps,dname)
%{ Calculate the material property target (MPT) curves for the 
% batch process. This function calculates both dimensionless
% curves and curves with dimensions.
%} 
    K_cross = 0;
    Q_cross = 0;
     
    v_mem = mpt_mass / rho_mem_SI; %[m3]
        
    N_BV = v_total / v_mem; %[-]
    
    K_bar = K_span * cin; %[-]
    
    r = cin/cout; %[-]
    
    % Calculate material property targets
    Q_ = (1/(1-eps)) * N_BV*(cin - cout) .* (1+K_span*cout)./(K_span*cout); %[mmol/cm3]
    
    Q = Q_./rho_mem; %[mmol/g]
   
    % Calculate dimensionless material property targets
    Q_nondim = nondim_bat(eps,N_BV,r,K_bar);
    
    % Predict cross over point
    [K_cross Q_cross] = crossover(eps,N_BV,r);
    
    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/batmpt_Q_',num2str(mpt_mass),'.csv'),Q);
    csvwrite(strcat(dname,'/batmpt_Qnondim_',num2str(mpt_mass),'.csv'),Q_nondim);
    csvwrite(strcat(dname,'/crossover',num2str(mpt_mass),'.csv'),[K_cross Q_cross]);
end