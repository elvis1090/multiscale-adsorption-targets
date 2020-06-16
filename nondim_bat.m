function [Q_bar_bat] = nondim_bat(eps,N_BV,rec,K_bar)
% This models the expression for dimensionless material property target 
% curves of the batch process
    Q_bar_bat = ((1/(1-eps))*N_BV*(1 - (1/rec)))*((rec + K_bar)./(K_bar));
end