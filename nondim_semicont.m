function [Q_bar_sb] = nondim_semicont(eps,tau,N_BV,K_bar)
% This models the expression for dimensionless material property target 
% curves of the semicontinuous process
    Q_bar_sb = (1/(1-eps))*(tau*N_BV-eps).*((1+K_bar)./K_bar);
end