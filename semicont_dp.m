function [dp] = semicont_dp(v_total,t_total,mu,m_mem,A_mem,del_P,N,rho_mem_SI)
% This function returns the membrane pore size required to operate at a
% given pressure drop (delta P), assuming that total mass of the membrane
% is contained within a single module (as opposed to multiple, parallel modules)

    v_mem = m_mem / rho_mem_SI;
    
    N_BV = v_total ./ v_mem;
    
    lp = v_mem / A_mem;
    
    dp = ((128*mu*lp.^2.*N_BV)./(pi*N*del_P*t_total)).^(1/4);           

end