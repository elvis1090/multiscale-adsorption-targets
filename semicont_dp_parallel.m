function [dp, n_mod, v_mem, lp] = semicont_dp_parallel(v_total,t_total,mu,m_mem,A_mem,del_P,N,rho_mem_SI,sp_thr)
% This function returns the membrane pore size required to operate at a
% given pressure drop (delta P), assuming that total mass of the membrane
% is split between mutliple paralle modules (as opposed to a single module
% with very high membrane thickness)

    v_mem = m_mem / rho_mem_SI;
    
    n_mod = ceil(v_mem/(sp_thr*A_mem));
    
    N_BV = v_total ./ v_mem;
    
    lp = v_mem ./ (n_mod*A_mem);
    
    dp = ((128*mu*lp.^2.*N_BV)./(pi*N*del_P*t_total)).^(1/4);
    
end