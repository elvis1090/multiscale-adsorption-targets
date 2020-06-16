function [K_ub,Q_ub] = semicont_capacity_ub_curve(v_total,t_total,t_regen,mu,m_mem,A_mem,del_P,N,rho_mem_SI,sp_thr,eps,rho_polym,N_A,v_bar,cin,rho_mem,dname)
%{ This function is used to calculate the upper bound curve for saturation 
% capacity of a membrane using the pore diameter model with parallel
% modules i.e. semicont_dp_parallel_sens
% NOTE: this function is superseded by the python/pyomo optimization script
% for maximizing saturation capacity. it is left in place as a dependancy
% for some functions such as the plotting code in the jupyter notebook.
%}
    fname_suff = del_P;
    
    % convert del_P to psi
    del_P = del_P * 6894.76;
    
    [dp_SI, n_mod, v_mem, lp] = semicont_dp_parallel(v_total,t_total,mu,m_mem,A_mem,del_P,N,rho_mem_SI,sp_thr); % dp in [m]
    
    dp = dp_SI * 1e2; % dp in [cm]

    Q_ub = semicont_maxcapacity(eps,rho_polym,dp,N_A,v_bar); % [mmol/g]

    [K_ub, N_BV] = semicont_mpt_inv(rho_mem_SI,m_mem,cin,rho_mem,Q_ub,v_total,t_total,t_regen,eps);

    writematrix([K_ub; Q_ub; N_BV; m_mem; dp; n_mod; v_mem; lp].',...
        strcat(dname,'/ub_curve_psi_',num2str(fname_suff),'.csv'));
end