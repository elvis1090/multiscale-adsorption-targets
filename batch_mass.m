function m_batch = batch_mass(v_total, cin, cout, Qmass, eps, K)
    % model function for the calculation of adsorbent mass when adsorption 
    % is done in batch mode given material
    % properites K, Q, epsilon and operating parameters v_total, cin, cout.

    % convert v_total to cm3
    v_total = v_total * 100^3;

    m_batch = (v_total*(cin - cout) / (Qmass*(1-eps))) * (1 + K*cout)/(K*cout); % [g]
end