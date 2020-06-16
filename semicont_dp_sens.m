function semicont_dp_sens(v_total,t_total,mu,m_mem,A_mem,del_P_span,N,rho_mem_SI,dname)
%{ This function implements the same function in semicont_porediameter, 
% but in vectorized form for a sensitivity analysis.
%}
    dp = semicont_dp(v_total,t_total,mu,m_mem,A_mem,del_P_span,N,rho_mem_SI);
        
    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/sbdpt_dp_',num2str(m_mem),'.csv'),dp);
end