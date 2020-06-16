# This script maximizes the upper bound on the saturation capacity

# load libraries
import numpy as np
from pyomo.environ import *
from pyo_q_ub_curve_init import get_capacity_upperbound_opt_params # easy data tranfer from the initialization code
from pyo_q_ub_curve_init import get_membrane_model_params # easy data tranfer from the initialization code

# maximize upper bound on saturation capacity
def capacity_ub_curve(delP_psi, v_total, water_years, K, init_vals, case):

    # flag for solution status
    solved = False

    # Unpack initial values
    initdp, initv_mem, initQmass, initn_mod\
    = [a for a in init_vals];

    # get or calculate membrane model parameters
    Np, mu, rho_mem, N_A, v_bar, eps, A_mem, sp_thr, cin = \
    get_membrane_model_params()
    delP = delP_psi * 6894.76; # [Pa]
    t_total = water_years*365*24*60*60; # [s]

    # calculate initial value for Qvol
    initQvol = initQmass * rho_mem

    # define optimization model
    model = ConcreteModel()

    # bounds on dp are used to constrain the feasible region to ensure convergence to a solution,
    # based on values from MATLAB simulations. dp bounds do not imply physical limitations / constraints
    # model.dp = Var(within=NonNegativeReals,bounds=(1e-10,2e-6),initialize=initdp) # [m]
    model.LOGdp = Var(within=Reals,bounds=(log(1e-10),log(2e-6)),initialize=log(initdp)) # [log(m)]
    model.v_mem = Var(within=NonNegativeReals,bounds=(0.001,120),initialize=initv_mem) # [m3]
    model.Qvol = Var(within=NonNegativeReals,bounds=(0.3,75),initialize=initQvol) # [mmol/cm3]
    # model.nmod = Var(within=NonNegativeReals,bounds=(0.1, nmod_ub), initialize=(nmod_ub-1)) # [-]
    model.l_mem = Var(within=NonNegativeReals,bounds=(sp_thr, sp_thr), initialize=(sp_thr)) # [m]

    print("intial values")
    print("LOGdp = ", value(model.LOGdp))
    print("v_mem = ", value(model.v_mem))
    print("Qvol = ", value(model.Qvol))
    print("l_mem = ", value(model.l_mem))
    # print("n_mod = ", value(model.nmod))

    # constraints
    def max_saturation_capacity_rule(m):
        return m.Qvol*100*exp(m.LOGdp)*(N_A*np.pi)**(1/3)*(3*v_bar)**(2/3) == 4 * eps * 4**(2/3)
    model.max_saturation_capacity_constraint = Constraint(rule=max_saturation_capacity_rule)

    def material_property_rule(m):
        return m.Qvol*(1-eps)*K*m.v_mem == (v_total-eps*m.v_mem)*(1+K*cin)
    model.material_property_constraint = Constraint(rule=material_property_rule)

    def pore_diameter_rule(m):
        # return (m.dp**4)*np.pi*Np*delP*t_total*(m.nmod**2)*(A_mem**2) == 128*mu*m.v_mem*v_total
        return 4*m.LOGdp == log(128*mu*v_total) - log(np.pi*Np*delP*t_total) + 2*log(m.l_mem) - log(m.v_mem)
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

    # extract optimum values
    # dp = value(model.dp); # [m]
    dp = exp(value(model.LOGdp)); # [m]
    v_mem = value(model.v_mem); # [m3]
    Qvol = value(model.Qvol); # [mmol/cm3]
    # n_mod = value(model.nmod); # [-]
    l_mem = value(model.l_mem); # [m]

    # convert saturation capacity to mass basis
    Qmass = Qvol / rho_mem; # [mmol/g]
    # calculate number of modules
    n_mod = v_mem / (l_mem*A_mem);
    # calculate number of bed volumes
    nbv = v_total / v_mem;

    # prevent saving nonsensical values if solution is not optimal
    if (not (results.solver.status == SolverStatus.ok and \
    results.solver.termination_condition == TerminationCondition.optimal)):
        dp = np.nan;
        v_mem = np.nan;
        Qmass = np.nan;
        Qvol = np.nan;
        n_mod = np.nan;
        l_mem = np.nan;
        nbv = np.nan;
        solved = False;
    else:
        print("optimum values")
        print("dp = ", dp)
        print("v_mem = ", v_mem)
        print("Qvol = ", Qvol)
        print("Qmass = ", Qmass)
        print("n_mod=",n_mod)
        print("l_mem=",n_mod)
        print("nbv=",nbv)
        solved = True;
    # prevent printing/saving nonsensical values if solution is not optimal

    # return dp, v_mem, Qvol, Q_mass, n_mod, l_mem, solved
    return dp, v_mem, Qvol, Qmass, n_mod, l_mem, nbv, solved;

# end: Function capacity_ub_curve()

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("Running NEW Optimization Problem")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

## Get code run parameters
v_total_all, water_years_all, cases, delP_psi_vec, n_el_mpt, notused1 = \
get_capacity_upperbound_opt_params()

# get membrane thickness
jnk1, jnk2, jnk3, jnk4, jnk5, jnk6, jnk7, sp_thr, jnk8 = \
get_membrane_model_params()

# vector to store K values across which sensitivity analyses was conducted
K_data = np.empty((n_el_mpt,len(delP_psi_vec)))

# vector to save optimum data
Q_ub_data = np.empty((n_el_mpt,len(delP_psi_vec)))
dp_data = np.empty((n_el_mpt,len(delP_psi_vec)))
n_mod_data = np.empty((n_el_mpt,len(delP_psi_vec)))
v_mem_data = np.empty((n_el_mpt,len(delP_psi_vec)))
l_mem_data = np.empty((n_el_mpt,len(delP_psi_vec)))
nbv_data = np.empty((n_el_mpt,len(delP_psi_vec)))

# vectors to store initial values, obtained from the solution of pyo_q_ub_curve_init.py
initQ_ub = np.zeros((n_el_mpt,len(delP_psi_vec)))
initdp = np.zeros((n_el_mpt,len(delP_psi_vec)))
# initn_mod = np.zeros((n_el_mpt,len(delP_psi_vec)))
initv_mem = np.zeros((n_el_mpt,len(delP_psi_vec)))

# default alternate initial values
# selected initial point is a solution from a succesful
# solve of the optimization problem
## set initial values [dp v_mem qmass nmod]
init_vals2 = [6.2e-9, 30.7, 3.8, 106]
init_vals3 = [1.68e-9, 14.94, 14.48, 519]
init_vals1 = [4.6e-9, 107.5, 5.3,108]

for case in cases:

    print("----------------------------------------------------------------------------------------------\
    ----------------------------------------------------------------------------------------------")
    print("CASE {0}".format(case))
    print("----------------------------------------------------------------------------------------------\
    ----------------------------------------------------------------------------------------------")

    dname = "case"+str(case)+"/max_cap_curve/sp_thr_"+str(format(sp_thr,'.0e'))
    #dname = "case"+str(case)+"/max_cap_curve/sp_thr_"+str(sp_thr)

    # read initial value data and vector of Ks for sensitivity analysis
    initQ_ub = np.genfromtxt(dname+"/pyo_cap_ub_Qs_init.csv",delimiter=',')
    initdp = np.genfromtxt(dname+"/pyo_cap_ub_dps_init.csv",delimiter=',')
    # initn_mod = np.genfromtxt(dname+"/pyo_cap_ub_nmods_init.csv",delimiter=',')
    initv_mem = np.genfromtxt(dname+"/pyo_cap_ub_v_mem_init.csv",delimiter=',')
    K_data = np.genfromtxt(dname+"/pyo_cap_ub_Ks.csv",delimiter=',')

    for i,delP_psi in enumerate(delP_psi_vec):

        for j,mpt_K in enumerate(K_data):

            print("**********************************************************************************************")
            print("K = {0} 1/M".format(mpt_K))
            print("**********************************************************************************************")

            ## set initial values [dp v_mem qmass nmod]
            # qmass will be converted to qvol in the optimization function
            init_vals = [initdp[j], initv_mem[j], initQ_ub[j], 9999999]

            if (init_vals[0] == 'nan' or np.isnan(init_vals[0])):
                if case == 1:
                    init_vals = init_vals1
                elif case == 2:
                    init_vals = init_vals2
                elif case == 3:
                    init_vals = init_vals3
                # end: switch over cases for initial values
            # end: initial value check and replacement

            dp, v_mem, Qvol, Qmass, n_mod, l_mem, nbv, solved = capacity_ub_curve(delP_psi, v_total_all[case-1], water_years_all[case-1], mpt_K, init_vals, case);

            if not solved: # try an alternate initial point.
                if case == 2:
                    # init_vals2 = [1.58e-9, 4.26, 15.47, 633]
                    dp, v_mem, Qvol, Qmass, n_mod, l_mem, nbv, solved = \
                    capacity_ub_curve(delP_psi, v_total_all[case-1], \
                    water_years_all[case-1], mpt_K, init_vals2, case);

                    if solved: # if successful, update alternate initial values
                        init_vals2 = [dp, v_mem, Qmass, n_mod]
                    # end: if successful, update alternate initial values

                elif case == 3:
                    # init_vals3 = [1.68e-9, 14.94, 14.48, 519]
                    dp, v_mem, Qvol, Qmass, n_mod, l_mem, nbv, solved = \
                    capacity_ub_curve(delP_psi, v_total_all[case-1], \
                    water_years_all[case-1], mpt_K, init_vals3, case);

                    if solved: # if successful, update alternate initial values
                        init_vals3 = [dp, v_mem, Qmass, n_mod]
                    # end: if successful, update alternate initial values

                elif case == 1:
                    # init_vals3 = [1.68e-9, 14.94, 14.48, 519]
                    dp, v_mem, Qvol, Qmass, n_mod, l_mem, nbv, solved = \
                    capacity_ub_curve(delP_psi, v_total_all[case-1], \
                    water_years_all[case-1], mpt_K, init_vals1, case);

                    if solved: # if successful, update alternate initial values
                        init_vals1 = [dp, v_mem, Qmass, n_mod]
                    # end: if successful, update alternate initial values

                # End: alternate starting initial values
            else:
                # if optimum solution is found, update alternate initial values
                # with nearest neighbor solution
                if case == 2:
                    init_vals2 = [dp, v_mem, Qmass, n_mod]
                elif case == 3:
                    init_vals2 = [dp, v_mem, Qmass, n_mod]
                elif case == 1:
                    init_vals1 = [dp, v_mem, Qmass, n_mod]
                # End nearest neighbor initialization update
            # End solution with alternate initial value

            # Group results into meaningful vectors (after trying alternate
            # initializations in case the first method fails)
            Q_ub_data[j,i] = Qmass; # [mmol/g]
            np.savetxt(dname+"/pyo_cap_ub_Qs.csv",Q_ub_data,delimiter=',');
            dp_data[j,i] = dp; # [m]
            np.savetxt(dname+"/pyo_cap_ub_dps.csv",dp_data,delimiter=',');
            n_mod_data[j,i] = n_mod; # [m]
            np.savetxt(dname+"/pyo_cap_ub_nmods.csv",n_mod_data,delimiter=',');
            v_mem_data[j,i] = v_mem; # [m]
            np.savetxt(dname+"/pyo_cap_ub_vmem.csv",v_mem_data,delimiter=',');
            l_mem_data[j,i] = l_mem; # [m]
            np.savetxt(dname+"/pyo_cap_ub_lmem.csv",l_mem_data,delimiter=',');
            nbv_data[j,i] = nbv; # [-]
            np.savetxt(dname+"/pyo_cap_ub_nbv.csv",nbv_data,delimiter=',');
        # end # loop over mpt_masses
    # end # loop over delP_psi

    # Summary file for easy analysis
    summary_array = np.array([K_data, Q_ub_data[:,0], dp_data[:,0], n_mod_data[:,0], v_mem_data[:,0], l_mem_data[:,0], nbv_data[:,0]])
    np.savetxt(dname+'/pyo_summary.csv',summary_array.T,delimiter=',',header="K(1/M),Q(mmol/g),dp(m),n_mod(-),v_mem(m3),l_mem(m),nbv(-)")

# end loop over case studies
