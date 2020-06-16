function nbv = cont_nbv(Kbar, Qbar, eps)
    % model function for the calculation of number of bed volumes when adsorption 
    % is done in semi-continuous mode given dimensionless material
    % properites Kbar, Qbar, epsilon
    
    nbv = Qbar * (1-eps) * ((Kbar)/(1+Kbar)) + eps;
end