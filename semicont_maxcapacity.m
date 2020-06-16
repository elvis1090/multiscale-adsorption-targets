function [Q] = semicont_maxcapacity(eps,rho_polym,dp,N_A,v_bar)
%{This function models the saturation capacity as a function of the material 
% structure, i.e. pore diameter - based on monolayer coverage and
% molar volume of the soute
%}
    Q = (eps/(1-eps))*(1/rho_polym)*(4./dp)*(1/(N_A*pi)^(1/3))*((4)/(3*v_bar))^(2/3); % [mmol/g]
end
