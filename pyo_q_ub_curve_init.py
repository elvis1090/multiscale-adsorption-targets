# This script generates initial values for the optimization problem to maximize
# the upper bound on the saturation capacity

from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import numpy as np
import sys

# get initial values for capacity upper bound optimization
def capacity_ub_curve_init(delP_psi, v_total, water_years, K, init_vals):

    # Unpack initial values
    initdp, initv_mem, initQvol = [a for a in init_vals];

    # get or calculate membrane model parameters
    Np, mu, rho_mem, N_A, v_bar, eps, A_mem, sp_thr, cin = \
    get_membrane_model_params()
    delP = delP_psi * 6894.76; # [Pa]
    t_total = water_years*365*24*60*60; # [s]

    model = ConcreteModel()

    # bounds on dp are used to constrain the feasible region to ensure convergence to a solution,
    # based on values from MATLAB simulations. dp bounds do not imply physical limitations / constraints
    model.dp = Var(within=NonNegativeReals,bounds=(1e-10,2e-6),initialize=initdp) # [m]
    model.v_mem = Var(within=NonNegativeReals,bounds=(1e-6,120),initialize=initv_mem) # [m3]
    model.Qvol = Var(within=NonNegativeReals,bounds=(0.05,75),initialize=initQvol) # [mmol/cm3]

    print("intial values")
    print("dp = ", value(model.dp))
    print("v_mem = ", value(model.v_mem))
    print("Qvol = ", value(model.Qvol))

    # NEW constraints
    def max_saturation_capacity_rule(m):
        return m.Qvol*100*m.dp*(N_A*np.pi)**(1/3)*(3*v_bar)**(2/3) == 4 * eps * 4**(2/3)
    model.max_saturation_capacity_constraint = Constraint(rule=max_saturation_capacity_rule)

    def material_property_rule(m):
        return m.Qvol*(1-eps)*K*m.v_mem == (v_total-eps*m.v_mem)*(1+K*cin)
    model.material_property_constraint = Constraint(rule=material_property_rule)

    def pore_diameter_rule(m):
        return (m.dp**4)*np.pi*Np*delP*t_total*m.v_mem == 128*mu*v_total*sp_thr**2
    model.pore_diameter_constraint = Constraint(rule=pore_diameter_rule)

    def obj_rule(m):
        return m.Qvol
    model.obj = Objective(rule=obj_rule,sense=maximize)

    solver = SolverFactory('ipopt',executable='/afs/crc.nd.edu/x86_64_linux/i/ipopt/3.12.8-hsl/bin/ipopt')
    solver.options['print_level'] = 0
    solver.options['linear_solver'] = 'ma27'
    solver.options['max_iter'] = 5000
    # solver.options['output_file'] = dname+str(K)+"ipopt.txt"
    # solver.options['file_print_level'] = 5
    results = solver.solve(model,tee=True)

    dp = value(model.dp); # [m]
    v_mem = value(model.v_mem); # [m3]
    Qvol = value(model.Qvol); # [mmol/cm3]

    # convert saturation capacity to mass basis
    Qmass = Qvol / rho_mem; # [mmol/g]

    # calculate number of modules
    n_mod = v_mem / (sp_thr*A_mem)

    # prevent saving nonsensical values if solution is not optimal
    if (not (results.solver.status == SolverStatus.ok and \
    results.solver.termination_condition == TerminationCondition.optimal)):
        dp = np.nan;
        v_mem = np.nan;
        Qmass = np.nan;
        Qvol = np.nan;
        n_mod = np.nan;
    else:
        print("optimum values")
        print("dp = ", dp)
        print("v_mem = ", v_mem)
        print("Qvol = ", Qvol)
        print("Qmass = ", Qmass)
        print("n_mod=",n_mod)
    # prevent printing/saving nonsensical values if solution is not optimal

    # return dp, v_mem, Qvol, Q_mass, n_mod
    return dp, v_mem, Qvol, Qmass, n_mod;

# end # Function capacity_ub_curve_init()

## read/set capacity upper bound optimization problem parameters
def get_capacity_upperbound_opt_params():

    # volume to be treated
    v_total_all = np.genfromtxt('v_total_all.csv',delimiter=',') # [m3]
    # treatment time
    water_years_all = np.genfromtxt('water_years_all.csv',delimiter=',') # [years]

    # case studies to run
    cases = [1,2,3]

    # pressure vector to calculate capacity upper bound curve on
    delP_psi_vec = [65]; # psi

    # number of points for upper bound curve
    n_el_mpt = 500;

    # vector of Ks for upper bound curve calculation
    cont_mpt_Ks = np.logspace(-1,8,num=n_el_mpt) # [1/M = l/mol = cm3/mmol]

    return v_total_all, water_years_all, cases, delP_psi_vec, n_el_mpt, cont_mpt_Ks
## END: read case study parameters

## set membrane model parameters
def get_membrane_model_params():

    Np = 1.6e14; # [1/m2] number of pores per m2 of membrane

    mu = 8.9e-4; # [Pa-s] viscosity of water at room temperature

    l_mem = 50e-6; # [m] membrane thickness
    # l_mem taken from value reported in the SI of DOI: 10.1021/acsami.7b04603

    # Area of 1 inch diameter of membrane
    A_mem_1 = (np.pi/4)*(2.54e-2)**2; # [m2]

    # Volume of 1 inch diameter of membrane
    vol_mem_1 = A_mem_1 * l_mem;  # [m3]

    # density of the membrane
    rho_mem = (17.6e-3)/(vol_mem_1*(1e2)**3); # [g/cm3]

    N_A = 6.022e20; # [1/mmol] Avogadro's number for 1 mmol of substance

    mw_Pb = 207.2e-3; # [g/mmol] Molecular mass of 1 mmol of lead i.e. target solute / contaminant

    rho_Pb = 11.34; # [g/cm3] Density of lead i.e. target solute / contaminant

    # molar volume of Pb [cm3/mmol]
    v_bar = mw_Pb/rho_Pb;

    # Note: hard coding void fraction / porosity based on a value Bill
    # mentioned during discussions on 03-Oct-2019. See Research meeting notes
    # on 04-Oct-2019 for more details.
    eps = 0.3; # [-]

    A_mem = ((np.pi/4)*(2.54e-2)**2)/(17.6e-6); # [m2] membrane area for 1 kg of membrane
    # Area (A_mem) calculated abased on information that 1" dia. of membrane weighs
    # 17.6 mg found in SI of doi: 10.1021/acs.langmuir.5b01605

    sp_thr = 1e-2 # [m] The maximum (threshold) value of membrane thickness for it to be formed into
    # a spiral wound module

    cin_ppb = 100; # [ppb] Lead concentration in Flint Michigan water
    # ref: 1) https://www.mlive.com/news/flint/2016/03/some_flint_water_test_sites_st.html
    # 2) doi: 10.1021/acsami.7b04603

    cin = cin_ppb * (1/207.2) / (1e3); #[mmol/cm3]

    return Np, mu, rho_mem, N_A, v_bar, eps, A_mem, sp_thr, cin
## END: set membrane model parameters

def main():

    print("initialization code")

    # default initial values. final initial values will be chosen lower in the code
    init_vals = [8.6E-10, 20.26, 19.7];

    ## Get problem parameters
    v_total_all, water_years_all, cases, delP_psi_vec, n_el_mpt, cont_mpt_Ks = \
    get_capacity_upperbound_opt_params()

    # get membrane thickness
    jnk1, jnk2, jnk3, jnk4, jnk5, jnk6, jnk7, sp_thr, jnk8 = \
    get_membrane_model_params()

    # vector to save data
    Q_ub_initdata = np.empty((n_el_mpt,len(delP_psi_vec)))
    dp_initdata = np.empty((n_el_mpt,len(delP_psi_vec)))
    n_mod_initdata = np.empty((n_el_mpt,len(delP_psi_vec)))
    v_mem_initdata = np.empty((n_el_mpt,len(delP_psi_vec)))

    for case in cases:

        print("----------------------------------------------------------------------------------------------\
        ----------------------------------------------------------------------------------------------")
        print("CASE {0}".format(case))
        print("----------------------------------------------------------------------------------------------\
        ----------------------------------------------------------------------------------------------")

        dname = "case"+str(case)+"/max_cap_curve/sp_thr_"+str(format(sp_thr,'.0e'))

        # vector of initial values, chosen based on value of membrane thickness (sp_thr)
        # (obtained from an early successful solve of the optimization problem)
        # [dp vmem qvol]
        if sp_thr == 1e-4:
            init_vals = [8.6E-10, 20.26, 19.7];
        elif sp_thr == 1e-2:
            init_vals = [5.3e-9, 59.35, 3.2]
        # end: switch initial values based on membrane thickness

        # save K data to files for plotting
        np.savetxt(dname+"/pyo_cap_ub_Ks.csv",cont_mpt_Ks,delimiter=',');

        for i,delP_psi in enumerate(delP_psi_vec):

            for j,mpt_K in enumerate(cont_mpt_Ks):

                print("**********************************************************************************************")
                print("K = {0} 1/M".format(mpt_K))
                print("**********************************************************************************************")

                dp, v_mem, Qvol, Qmass, n_mod = capacity_ub_curve_init(delP_psi, v_total_all[case-1], water_years_all[case-1], mpt_K, init_vals);

                # Group results into meaningful vectors
                Q_ub_initdata[j,i] = Qmass; # [mmol/g]
                np.savetxt(dname+"/pyo_cap_ub_Qs_init.csv",Q_ub_initdata,delimiter=',');
                dp_initdata[j,i] = dp; # [m]
                np.savetxt(dname+"/pyo_cap_ub_dps_init.csv",dp_initdata,delimiter=',');
                n_mod_initdata[j,i] = n_mod; # [m]
                np.savetxt(dname+"/pyo_cap_ub_nmods_init.csv",n_mod_initdata,delimiter=',');
                v_mem_initdata[j,i] = v_mem; # [m]
                np.savetxt(dname+"/pyo_cap_ub_v_mem_init.csv",v_mem_initdata,delimiter=',');

                # replace initial values only if previous problem was solved to
                # optimality
                if not(np.isnan(dp)):
                    init_vals = [dp, v_mem, Qvol];
                # end: replace initial values only if previous problem was solved to
                # optimality

            # end # loop over mpt_masses

        # end # loop over delP_psi

        # Summary file for easy analysis
        summary_array = np.array([cont_mpt_Ks, Q_ub_initdata[:,0], dp_initdata[:,0], n_mod_initdata[:,0], v_mem_initdata[:,0]])
        np.savetxt(dname+'/pyo_summary_init.csv',summary_array.T,delimiter=',',header="K_init(1/M),Q_init(mmol/g),dp_init(m),n_mod_init(-),v_mem_init(m3)")

    # end loop over case studies

# end Main

# run main code
if __name__ == "__main__":
    main();
# END: run main code
