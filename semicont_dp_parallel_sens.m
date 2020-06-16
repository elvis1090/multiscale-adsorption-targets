function semicont_dp_parallel_sens(v_total,t_total,mu,m_mem,A_mem,del_P_span,N,rho_mem_SI,sp_thr,dname)
%{ Pore diameter targets for multiple membrane modules in parallel wherein
% number of modules is determined by the function such that the thickness 
% of the membrane in each module is <= 1mm. Membrane mass, and hence 
% thickness is uniform across all modules. Function is a modified form of 
% semicont_dptargets. This function is used for sensitivity analyses
%}

    [dp, n_mod] = semicont_dp_parallel(v_total,t_total,mu,m_mem,A_mem,...
        del_P_span,N,rho_mem_SI,sp_thr);

    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/sbdpt_dp_parallel',num2str(m_mem),'.csv'),dp);
    csvwrite(strcat(dname,'/sbdpt_nmod_parallel',num2str(m_mem),'.csv'),n_mod);

end