function calc_dmlss_props(dname)
    % This function is used to calculate the dimensionless properties K and Q
    % for analysis (or plotting) with our dimensionless plots. This
    % function also calculates and saves into the same file the number of 
    % bed volumes that can be treated for the candidate material in the
    % arsenic removal case study by a semi-continuous adsorptive process
    
    disp(' ');
    disp('Calculating dimensionless properties for arsenic case study');

    % read material property data
    arsenic_ads_dat = readmatrix('./arsenic_adsorbents.csv',...
        'Range','B:C'); % [l/g_arsenic, g_arsenic/g_adsorbent]
    
    % get number of data points
    n = length(arsenic_ads_dat);
    
    % initialize an empty vector to store 
    nbv_all = zeros(n,2);
    
    % calculation parameters
    % NOTE: cin_gl and cout_gl should correspond to the same recovery
    % ratio in the last element of the rec_span array (and thereby associated parameters) 
    % in adsorption_analysis
    % for example: cin_gl = 140e-6, cout_gl = 10e-6 ==> r = cin_gl/cout_gl
    % = 14. last element of rec_span = 14
    cin_gl = 140e-6; % [g_arsenic/l], Sarkar et al. 10.1016/j.watres.2010.07.072
    cout_gl = 10e-6; % [g_arsenic/l], US-EPA https://www.epa.gov/ground-water-and-drinking-water/national-primary-drinking-water-regulations#seven.
    rho_ads = 700; % [g_adsorbent/l_adsorbent], assumed 0.7 g/cm3
    eps = 0.3; % assumed, average value for porous adsorbents
    r = cin_gl/cout_gl
    
    % convert Q from mass to volume basis
    Q_arsenic_mass = arsenic_ads_dat(:,2) * rho_ads;    
    
    % calculate
    Kbar_arsenic = arsenic_ads_dat(:,1)*cin_gl;
    Qbar_arsenic = Q_arsenic_mass/cin_gl;
    
    % save data to meaningful variables for easy storage
    K_arsenic = arsenic_ads_dat(:,1);
    Q_arsenic = arsenic_ads_dat(:,2);
    
    
    % calculate number of bed volumes that can be treated
    for i = 1:n
        nbv_all(i,1) = cont_nbv(Kbar_arsenic(i), Qbar_arsenic(i), eps);
        nbv_all(i,2) = batch_nbv(Kbar_arsenic(i), Qbar_arsenic(i), eps, r);
    end 
    
    % write data to file to plot in python
    writematrix([K_arsenic, Q_arsenic, Kbar_arsenic, Qbar_arsenic, nbv_all],strcat(dname,'/arsenic_adsorbents_dimensionless.csv'));
    
end