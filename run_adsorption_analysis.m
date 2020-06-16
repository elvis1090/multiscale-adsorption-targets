% This code is used to perform analysis of the adsorptive membrane process
% for three different cases: 
% Case 1: 60 m3 in 2 years. This is the baseline case where 60m3 of water 
% is the average potable water consumption of a small household in 2 years.
% Case 2: 15 m3 in 6 months (0.25*2 years). Regeneration is carried out 
% every 6 months.
% Case 3: (Not used in paper) 15m3 in 2 years. The idea here is that we will have parallel
% trains installed, each treating about a quarter of the two-year
% consumption of the household. This design enables a lower interstitial
% velocity through the membrane, allowing the use of smaller pore sizes
% which could raise the upper bound on saturation capacity and possibly
% enable lower requirement of membrane masses for a process
% Case 4: Sensitivity analysis over dimensionless material property curves
% for relative behavior of batch and semi-continuous processes and arsenic
% case study and dimensionless properties of existing materials for arsenic 
% removal case study
% Case 5: Analyze over-/under-predictions due to incorrect scale up.
% Case 6: Supporting information part 1 - breakthrough time calculations, 
% figures for derivation for competitive behavior
% between batch and semi-continuous process configurations for adsorption
% Case 7: Supporting information part 2 - Nonlinear regression for Langmuir
% isotherm
% Case 8: Supporting information part 3 - Reference dimensionless plots

% clear outputs and memory
clear all; close all; clc;

% vector of cases for iteration
% choose which case needs to be run
% cases = [1,2,3,4,5,6,7,8] runs all analysis including those in the
% supporting information
cases = [1,2,3,4,5,6,7,8]; % corresponding to descriptions in header

% vector of water consumption corresponding to cases
v_total_all = [60 15 15]; %[m3]
% save to file to synchronize with python scripts (capacity upper bound)
% curve calculation
writematrix(v_total_all,'v_total_all.csv');

% vector of time in years for water consumption corresponding to above
% cases
water_years_all = [2 0.25*2 2]; %[years]
% save to file to synchronize with python scripts (capacity upper bound)
% curve calculation
writematrix(water_years_all,'water_years_all.csv');

% vector of regeneration times corresponding to above cases 
% for all foreseeable applications, this is the same as consumption time
regen_years_all = water_years_all; %[years]

% vector of module masses for breakthrough time calculation corresponding
% to above cases
module_mass_all = [46 50 50]; %[kg]

% number of cases
n_cases = length(cases);

% iterate over cases and perform analysis
for i = 1:n_cases
    % define name of case study
    casename = strcat('case',num2str(cases(i)));
    disp('------------------------------------------------------');
    disp(['                   ',casename,'                     ']);
    disp('------------------------------------------------------');
    % directory to save results
    dname = casename;
    if (strcmp(dname,'case4'))
        v_total = v_total_all(2);
        t_total = water_years_all(2)*365*24*60*60;
        t_regen = regen_years_all(2)*365*24*60*60;
        module_mass = module_mass_all(1);        
    elseif(strcmp(dname,'case5') || strcmp(dname,'case6')||...
            strcmp(dname,'case7') || strcmp(dname,'case8'))
        % note that analyzing over-/under-predictions, 
        % nonlinear regression of Langmuir isotherm and 
        % dimensionless sensitivity analysis are independent of the
        % specific values of v_total, t_total, etc. chosen and is hence
        % clubbed with breakthrough time calculations for POU_24^C
        v_total = v_total_all(1);
        t_total = water_years_all(1)*365*24*60*60;
        t_regen = regen_years_all(1)*365*24*60*60;
        module_mass = module_mass_all(1); 
    else
        v_total = v_total_all(cases(i));
        t_total = water_years_all(cases(i))*365*24*60*60; %[s]
        t_regen = regen_years_all(cases(i))*365*24*60*60; %[s]
        module_mass = module_mass_all(cases(i)); %[kg]
    end
    adsorption_analysis(v_total,t_total,t_regen,module_mass,dname);
end

