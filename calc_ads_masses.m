function calc_ads_masses(v_total, cin, eps, rho_mem, cout, dname)
    % This function iterates over properties of existing adsorbents and
    % calculates the mass of membrane required for a given process defined by
    % v_total

    % read data and group into meaningful variables
    QK_dat = readmatrix('./sorbent_prop.csv','Range','B:C');
    Q_dat = QK_dat(:,1); %[mmol/g]
    K_dat = QK_dat(:,2); %[cm3/mmol]

    % get number of data points
    len = length(QK_dat);

    % initialize empty vectors to store masses for batch and semi-continuous
    % adsorption processes
    masses_batch = zeros(len,1);
    masses_cont = zeros(len,1);

    % iterate over adsorbents and calculate masses
    for i = 1:len
        masses_batch(i) = batch_mass(v_total, cin, cout, Q_dat(i), eps, K_dat(i)); %[g]
        masses_cont(i) = cont_mass(v_total, cin, Q_dat(i), eps, K_dat(i), rho_mem); %[g]
    end
    
    % convert masses to kg
    masses_batch = masses_batch/1000;
    masses_cont = masses_cont/1000;
    
    % save results to .csv file
    writematrix([K_dat, Q_dat, masses_batch, masses_cont],strcat(dname,'/adsorbent_masses.csv'));
end