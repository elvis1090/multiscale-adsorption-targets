function [dc_dt] = semicont_ode(t,c,v_total,t_total,eps,A_mem,nu,K,Q,cin,n_el,del_z)
%{
 A function that models the behaviour of fixed bed adsorption based on equilibrium wave theory.
 Ref:
 https://www.annualreviews.org/doi/full/10.1146/annurev-chembioeng-061312-103318
 OR doi: 10.1146/annurev-chembioeng-061312-103318
 Eq. 5 from the paper has been modeled by performing a finite element
 discretization of the space (z) variable using the backward difference
 method. The formula implemented below appears flipped only due a the sign convention
 of the conservation equations.

Inputs:
V: interstitial velocity [m/s]
nu: phase ratio [dimensionless]
K: binding affinity [cm3/mmol]
Q: saturation capacity [mmol/cm3]
cin: concentration of fluid at the moment adsorption starts [mmol/cm3]
n_el: number of elements in finite element discretization [dimensionless]
del_z: differential length [m]

Outputs:
dc_dt: vector for numerical integration
 %}
    
    dc_dt = zeros(n_el,1);
    
    V = (v_total)/(t_total*eps*A_mem);

    for i = 1:n_el
        if (i==1)
            dc_dt(i) = (V/(1+nu*K*Q/(1+K*c(i))^2))*(cin - c(i))/del_z;
            % modeling as above i.e. cin - c(i) for the first element and
            % c(i-1) - c(i) (below) propagates the plug / shock of change
            % of concentrations when adsorption just starts
        else
            dc_dt(i) = (V/(1+nu*K*Q/(1+K*c(i))^2))*(c(i-1) - c(i))/del_z;
        end
    end
    

end