function [K,N_BV] = semicont_mpt_inv(rho_mem_SI,mpt_mass,cin,rho_mem,Q_span,v_total,t_total,t_regen,eps)
%{This function returns K from the model for the mpt curves given (fixed)
% mpt_mass and Q. 
%}   
    v_mem = mpt_mass / rho_mem_SI; %[m3]
    
    N_BV = v_total ./ v_mem; %[-]
    
    tau = t_regen/t_total; %[-]
    
    logK = -log(Q_span*rho_mem*(1-eps)./(tau*N_BV-eps) - cin);
    
    % get rid of comlex numbers that arise due to taking logarithm of a
    % negative number
    logK(imag(logK)~=0) = -999; 

    K = exp(logK);
end