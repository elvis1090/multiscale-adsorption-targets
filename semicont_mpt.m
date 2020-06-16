function semicont_mpt(rho_mem_SI,mpt_mass,cin,rho_mem,K_span,v_total,t_total,t_regen,eps,dname)
%{ Calculate the material property target (MPT) curves for the 
% semicontinuous process. This function calculates both dimensionless
% curves and curves with dimensions.
%} 
    K = K_span;
    
    v_mem = mpt_mass / rho_mem_SI; %[kg/m3]
    
    N_BV = v_total / v_mem; %[-]
    
    tau = t_regen / t_total;
    
    K_bar = K_span*cin; %[-]
   
    Q_ = (1/(1-eps))*(tau*N_BV - eps).*((1+K*cin)./K);
    
    % convert Q_ to mmol/g units
    Q = Q_./rho_mem; %[mmol/g]
    
    % calculate dimensionless material property curves
    Q_nondim = nondim_semicont(eps,tau,N_BV,K_bar);
    
    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/sbmpt_Q_',num2str(mpt_mass),'.csv'),Q);
    csvwrite(strcat(dname,'/sbmpt_Qnondim_',num2str(mpt_mass),'.csv'),Q_nondim);
end