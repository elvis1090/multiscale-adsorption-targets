function nbv = batch_nbv(Kbar, Qbar, eps, r)
    % model function for the calculation of number of bed volumes when adsorption 
    % is done in batch mode given dimensionless material
    % properites Kbar, Qbar, epsilon, r
    
    nbv = (1-eps)*Qbar*Kbar/((1-1/r)*(r+Kbar));
end