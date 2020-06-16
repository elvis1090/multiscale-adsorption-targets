function [xK xQ] = crossover(eps, N_BV, r)
%{ Calculate the crossover point seen in Type 2 relative 
% behavior between batch and semi-continuous processes.
%}
    xK = (eps - N_BV*(2-r))./(N_BV/r - eps);
    xQ = nondim_bat(eps,N_BV,r,xK);
end