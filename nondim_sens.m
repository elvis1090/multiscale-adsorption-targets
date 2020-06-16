function nondim_sens(eps,N_BV,rec,K_bar,tau,dname)
% for sensitivity analysis of the nondimensional material property curves
% checking for sensitivity to eps and r

    % batch 
    Q_nd_bat = nondim_bat(eps,N_BV,rec,K_bar);
    
    % semi-continuous
    Q_nd_semi = nondim_semicont(eps,tau,N_BV,K_bar);
    
    [K_cross Q_cross] = crossover(eps,N_BV,rec);
    
    % save outputs to a .csv file (to plot in python)
    csvwrite(strcat(dname,'/Q_nd_bat_eps',num2str(eps),...
        '_r',num2str(rec,'%.1f'),'nbv_',num2str(N_BV,'%.2f'),'.csv'),Q_nd_bat.');
    csvwrite(strcat(dname,'/Q_nd_semi_eps',num2str(eps),...
        '_r',num2str(rec,'%.1f'),'nbv_',num2str(N_BV,'%.2f'),'.csv'),Q_nd_semi.');
    csvwrite(strcat(dname,'/crossover_eps',num2str(eps),...
        '_r',num2str(rec,'%.1f'),'nbv_',num2str(N_BV,'%.2f'),'.csv'),[K_cross Q_cross]);
end