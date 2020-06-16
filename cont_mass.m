function m_cont= cont_mass(v_total, cin, Qmass, eps, K, rho_mem)
    % model function for the calculation of adsorbent mass when adsorption 
    % is done in semi-continuous mode given material
    % properites K, Q, epsilon, rho_mem and operating parameters v_total, cin

    % convert v_total to cm3
    v_total = v_total * 100^3;

    m_cont = v_total / (((1-eps)*Qmass*K/(1+K*cin)) + eps/rho_mem); % [g]

end