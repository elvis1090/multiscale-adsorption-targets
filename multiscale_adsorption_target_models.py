'''
This script defines functions needed to caculate material property targets and
capacity upper bound curve for membrane adsorbents.

Elvis A. Eugene, William A. Phillip, Alexander W. Dowling
University of Notre Dame
'''

'''
Pseudocode / workflow for capacity upper bound curve calculation.
See Qmax_modified_methods.ipynb, Section 3: Updated method 2
'''

# load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import itertools
from scipy.optimize import fsolve
# from pyo_q_ub_curve_init import get_membrane_model_params # membrane parameters that were used in MATLAB code

## ---------------------------------------------------------------------------------------------- ##
##                                       helper functions                                         ##
## ---------------------------------------------------------------------------------------------- ##

# row major legend
# https://stackoverflow.com/questions/10101141/matplotlib-legend-add-items-across-columns-instead-of-down
def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])
# END flip()

# https://stackoverflow.com/questions/27579898/how-to-create-square-log-log-plots-in-matplotlib
def set_aspect_ratio_log(plot, aspect_ratio):
    '''
    inputs
    plot: matplotlib axis object
    aspect_ratio: final aspect ratio of plot, float scalar

    returns
    plot resized based on aspect_ratio
    '''
    x_min, x_max = plot.get_xlim()
    y_min, y_max = plot.get_ylim()
    return plot.set_aspect(aspect_ratio * ((np.log10(x_max / x_min)) / (np.log10(y_max / y_min))))
# END set_aspect_ratio_log()

## ------------------------------------ END helper functions ------------------------------------ ##

## ---------------------------------------------------------------------------------------------- ##
##                                       model functions                                          ##
## ---------------------------------------------------------------------------------------------- ##
def QB_target(K, NBV, r, epsilon=0.3):
    ''' Dimensionless adsorption capacity material property target for batch process

    Arguments:
        K: dimensionless binding affinity [float or vector]
        NBV: number of bed volumes treated before regeneration [float]
        r: rejection ratio [need to check precise name] [float]
    '''

    return (1/(1-epsilon))*NBV*(1-1/r)*(r+K)/K
# END QB_target()

def QC_target(K, NBV, epsilon=0.3):
    ''' Dimensionless adsorption capacity material property target for semi-continous process

    Arguments:
        K: dimensionless binding affinity [float or vector]
        NBV: number of bed volumes treated before regeneration, dimensionless [float]
    '''

    return (1/(1-epsilon))*(NBV - epsilon)*(1+K)/K
# END QC_target()

def pore_size_limit(N_BV,t_total,mu,l_mem,N_p,delP):
    '''
    inputs
    N_BV: vector of bed volumes [-], scalar or 1D numpy array
    t_total: regeneration time [s], scalar
    mu: fluid viscosity [Pa s], scalar
    l_mem: membrane thickness [m], scalar
    N_p: number of pores per unit area of membrane [1/m^2], scalar
    delP: available pressure drop [Pa], scalar

    returns
    dp: pore size [m], scalar or 1D numpy based on N_BV
    '''

    dp = ((N_BV/t_total)*((128*mu*l_mem**2)/(N_p*np.pi*delP)))**(1/4)

    return dp
# END pore_size_limit()

def saturation_limit(eps,dp,N_A,v_s_bar):
    '''
    inputs
    eps: porosity [-], scalar
    dp: pore size [m], scalar or 1D numpy array
    N_A: Avogadro's number [1/mmol], scalar
    v_s_bar: molar volume of solute [cm^3/mmol], scalar

    returns
    qmax: maximum saturation capacity of membrane [mmol/cm^3], scalar or 1D numpy based on N_BV
    '''

    qmax = (4*eps/((1-eps)*dp*1e2))*(1/(np.pi*N_A)**(1/3))*(4/(3*v_s_bar))**(2/3)

    return qmax
# END saturation_limit()

def calc_K(N_BV,Q,eps,cin):
    '''
    inputs:
    N_BV: number of bed volumes [-], 1D numpy array
    Q: Saturation cacpcity [mmol/cm^3], 1D numpy array
    eps: porosity [-], scalar
    cin: inlet concentration [mmol/cm^3], scalar

    returns:
    K: binding affinity [cm^3/mmol], 1D numpy array
    '''

    K = (N_BV-eps) / (Q*(1-eps)-cin*(N_BV-eps))

    return K
# END calc_K()

def bed_volume_limit(N_BV,eps,cin,t_total,N_p,delP,mu,l_mem,N_A,v_s_bar):
    '''
    N_BV: bed volumes [-], scalar or 1D numpy array
    eps: porosity [-], scalar
    cin: inlet concentration [mmol/cm^3], scalar
    t_total: regeneration time [s], scalar
    N_p: number of pores per unit area of membrane [1/m^2], scalar
    delP: available pressure drop [Pa], scalar
    mu: fluid viscosity [Pa s], scalar
    l_mem: membrane thickness [m], scalar
    N_A: Avogadro's number [1/mmol], scalar
    v_s_bar: molar volume of solute [cm^3/mmol], scalar

    returns
    res: residual of bed volume limit equation to solve with a nonlinear solver
    '''

    # split the equation into 4 terms from left to right for convenience
    term1 = 0.04*eps/cin
    term2 = ((t_total/N_BV)*((N_p*np.pi*delP)/(128*mu*l_mem**2)))**(1/4)
    term3 = (1/(np.pi*N_A))**(1/3)
    term4 = (4/(3*v_s_bar))**(2/3)

    res = term1 * term2 * term3 * term4 + eps - N_BV

    return res
# END bed_volume_limit()

def calc_semicont_masses(v_total, rho_ads, eps, Qvol, K, cin):
    '''
    Calculates the amount of adsorbent needed for a given separation

    inputs
    v_total: volume of water to be treated, float scalar [cm3]
    rho_ads: density of membrane, float scalar [g/cm3_mem]
    eps: (bed) porosity / void fraction, float scalar [unitless or cm3_pores/cm3_adsorbent bed]
    Qvol: saturation capacity, float numpy array [mmol/cm3_ads]
    K: binding affinity [cm3_solute/mmol]
    cin: concentration of solute in water [mmol/cm3_solute]

    outputs
    m_ads: mass of adsorbent needed for separation, float numpy array [g]
    '''

    m_ads = (v_total*rho_ads) / (((1-eps)*Qvol*K)/(1+K*cin) + eps) # [g]

    return m_ads
# end calc_semicont_masses()

def calc_batch_masses(v_total, rho_ads, eps, Qvol, K, cin, cout):
    '''
    Calculates the amount of adsorbent needed for a given separation

    inputs
    v_total: volume of water to be treated, float scalar [cm3]
    rho_ads: density of membrane, float scalar [g/cm3_mem]
    eps: (bed) porosity / void fraction, float scalar [unitless or cm3_pores/cm3_adsorbent bed]
    Qvol: saturation capacity, float numpy array [mmol/cm3_ads]
    K: binding affinity [cm3_solute/mmol]
    cin: concentration of solute in water [mmol/cm3_solute]

    outputs
    m_ads: mass of adsorbent needed for separation, float numpy array [g]
    '''

    m_batch = (v_total*rho_ads*(cin - cout) / (Qvol*(1-eps))) * (1 + K*cout)/(K*cout) # [g]

    return m_batch
# end calc_semicont_masses()


# ----------- packed bed models ---------------#
def bed_volumes(D_AB,t_total,Dp,lub_bar=0.01):
    '''
    Returns the number of bed volumes given diffusivity, particle size, and time for treatment
    Useful to calculate range of n_bv given Dp_min and Dp_max

    inputs
    D_AB : solute diffusivity [cm2/s]
    t_total : time available for treatment [s]
    d_particle : diameter of particle [cm]
    lub_bar : fractional length of unused bed (LUB/lb), default 0.01 [-]

    outputs
    n_bv : number of bed volumes treated [unitless]
    '''

    n_bv = lub_bar * D_AB * t_total / (Dp**2)

    return n_bv
# bed_volumes()

def bed_dimensions(v_total,n_bv,lb_db_ratio=3):
    '''
    calculates the length of the packed bed based on operating parameters and assumptions

    inputs
    v_total : volume of water to be treated [cm3]
    n_bv : number of bed volumes to be treated [unitless*]
    lb_db_ratio : ratio of length of bed to diameter of bed, default = 3

    outputs
    d_b : diameter of packed bed [cm]
    l_b : length of packed bed [cm]
    n_bed : number of adsorbent beds in paralle
    '''

    n_bed = np.ones(n_bv.shape) # we start with a single bed

    d_b = (4 * v_total / (lb_db_ratio * np.pi * n_bv))**(1/3) # [cm]

    l_b = lb_db_ratio * d_b # [cm]

    # if bed diameter exceeds 4.5 m, consider multiple beds in parallel by fixing d_b and l_b
    if d_b.max() >= 450:

        # find all indices where diameter is greater than 4.5 m  to adjust dimensions using
        # parallel beds
        over_limit_idx = d_b>450
#         print('over_limit_idx=',over_limit_idx)

        d_b_max = 450 # [cm] max bed diameter is 4.5 m from literature

        l_b_max = lb_db_ratio * d_b_max # [cm]

        v_bed_max = np.pi * d_b_max**2 * l_b_max / 4 # [cm3]

        v_ads = v_total/n_bv[over_limit_idx]

        n_bed[over_limit_idx] = v_ads / v_bed_max

        d_b[over_limit_idx] = (4 * v_total / (n_bed[over_limit_idx] * lb_db_ratio * np.pi * n_bv[over_limit_idx]))**(1/3)

        l_b[over_limit_idx] = lb_db_ratio * d_b[over_limit_idx]
    # END if bed diameter exceeds 4.5 m, consider multiple beds in parallel by fixing d_b and l_b

    return d_b, l_b, n_bed
# bed_dimensions()

def particle_pore_diameter(eps_part, n_ps, d_particle):
    '''
    calculates the pore diameter of a particle in the packed bed given operating conditions and other bed design criteria

    inputs
    eps_part : particle porosity [unitless*]
    n_ps : number of pores per spherical adsorbent particle [unitless]
    d_particle : diameter of the particle [cm]

    outputs
    d_pore : diameter of pore [cm]
    '''

    d_pore = (2*eps_part*d_particle**2 / (3*n_ps))**(1/2) # [cm]

    return d_pore
# particle_pore_diameter()

def saturation_capacity(eps_part, d_pore, n_a, v_s_bar):
    '''
    calculates the saturation capacity of the adsorbent based on particle and solute characteristics

    inputs
    eps_part : particle porosity [unitless*]
    d_pore : diameter of pore [cm]
    n_a : Avogadro's number [1/mmol]
    v_s_bar : molar volume of solute [cm3/mmol]

    outputs
    q_max : saturation capacity of the adsorbent [mmol/cm3(particle)]
    '''

    q_max = (4*eps_part/d_pore) * (1/(n_a*np.pi))**(1/3) * (4/(3*v_s_bar))**(2/3) # [mmol/cm3(particle)]

    return q_max
# saturation_capacity()

def binding_affinity(n_bv, eps_bed, Q, c_in):
    '''
    calculates the binding affinity corresponding to the saturation capacity based
    on the Langmuir isotherm and conservation equation for the semi-continuous adsorption process
    (duplicate of calc K?)

    inputs
    n_bv : number of bed volumes treated [unitless]
    eps_bed : bed void fraction [unitless]
    Q : saturation capacity [mmol/cm3(particle)]
    c_in : inlet concentration of solute [mmol/cm3(solute)]

    outputs
    K : binding affinity [cm3(solute)/mmol]
    '''

    K = (n_bv - eps_bed) / ((1-eps_bed)*Q - c_in*(n_bv - eps_bed))

    return K
# binding_affinity()

def particle_diameter(D_AB, t_total, n_bv,lub_bar=0.01):
    '''
    calculates diameter of particle given diffusivity of solute and number of bed volumes

    inputs
    D_AB : solute diffusivity [cm2/s]
    t_total : time available for treatment [s]
    n_bv : number of bed volumes treated [unitless]
    lub_bar : fractional length of unused bed (LUB/lb), default 0.01 [-]

    outputs
    d_particle : diameter of particle [cm]
    '''

    d_particle = (lub_bar*D_AB*t_total/(n_bv))**(1/2)

    return d_particle
# particle_diameter()

def pressure_drop(l_b, mu, eps_bed, n_bv, d_particle, t_total):
    '''
    calculates the pressure drop across the bed using the Kozeny Carman equation

    inputs
    l_b : length of packed bed [cm]
    mu : viscosity of solvent [Pa s]
    eps_bed : porosity or void fraction of the bed [unitless*]
    n_bv : number of bed volumes to be treated [unitless*]
    d_particle : diameter of the particle [cm]
    t_total : time available for treatment [s]

    outputs
    delP: pressure drop in packed bed [Pa]
    '''

    delP = l_b**2 * 150 * mu * (1-eps_bed)**2 * n_bv / (d_particle**2 * eps_bed**3 * t_total) # [Pa]

    return delP
# pressure_drop()

def peclet_number(d_particle,v_total,n_bed,D_AB,d_b,t_total):
    '''
    calculates the peclet number given operating and design parameters of the packed bed

    inputs
    d_particle : diameter of the particle [cm]
    v_total : volume to be treated [cm3]
    n_bed : number of beds in parallel [unitless]
    D_AB : solute diffusivity [cm2/s]
    d_b : diameter of packed bed [cm]
    t_total : time available for treatment [s]

    outputs
    n_pe : peclet number [dimensionless]
    '''

    n_pe = d_particle*4*v_total/(n_bed*D_AB*np.pi*d_b**2*t_total)

    return n_pe
# peclet_number()

def nbv_semicont(K,Q,eps_bed,cin):
    '''
    calculates the no. of bed volumes that can be treated given Langmuir isotherm
    constants and operating parameters

    inputs
    K : binding affinity [cm3/mmol]
    Q : saturation capacity [mmol/cm3]
    eps_bed : bed void fraction [cm3(voids)/cm3(bed)]
    cin : inlet concentration of solute [mmol/cm3]

    outputs
    nbv : number of bed volumes that can be treated in the semicontinuous mode [-]
    '''

    nbv = ((1-eps_bed) * Q * K / (1 + K*cin)) + eps_bed

    return nbv
# nbv_semicont()

def lub_pore_diameter(Q, eps_part, n_a, v_s_bar, d_pore_min_max=[None,None]):
    '''
    calculates the pore diameter of the adsorbent particle based on saturation capacity

    inputs
    eps_part : particle porosity [unitless*]
    n_a : Avogadro's number [1/mmol]
    v_s_bar : molar volume of solute [cm3/mmol]
    Q : saturation capacity of the adsorbent [mmol/cm3(particle)]
    d_pore_min_max : lower and upper limits for pore diameter [cm], list floats, default [None(min),None(max)]

    outputs
    d_pore : diameter of pore [cm]
    '''
    d_pore = (4*eps_part/Q) * (1/(n_a*np.pi))**(1/3) * (4/(3*v_s_bar))**(2/3) # [cm]

    # bound particle pore diameters if provided
    if d_pore_min_max[0]:
        d_pore[d_pore<d_pore_min_max[0]] = d_pore_min_max[0]
        d_pore[d_pore>d_pore_min_max[1]] = d_pore_min_max[1]
    # End bounding particle pore diameter

    return d_pore
# lub_pore_diameter()

def lub_particle_diameter(eps_part, n_ps, d_pore):
    '''
    calculates the pore diameter of a particle in the packed bed given operating conditions and other bed design criteria

    inputs
    eps_part : particle porosity [unitless*]
    n_ps : number of pores per spherical adsorbent particle [unitless]
    d_pore : diameter of pore [cm]

    outputs
    d_particle : diameter of the particle [cm]
    '''

    d_particle = d_pore*((3*n_ps) / (2*eps_part))**(1/2) # [cm]

    return d_particle
# lub_particle_diameter()

def lub_lub_bar(d_part,nbv,D_AB,t_total):
    '''
    calcualtes lub_bar given operating and design conditions

    inputs
    D_AB : solute diffusivity [cm2/s]
    t_total : time available for treatment [s]
    d_part : diameter of particle [cm]
    n_bv : number of bed volumes treated [unitless]

    outputs
    lub_bar : fractional length of unused bed (LUB/lb), default 0.01 [-]
    '''
    lub_bar = nbv * d_part**2 / (D_AB * t_total)

    return lub_bar
# lub_lub_bar()

def lub_bed_dimensions(v_total,nbv,dbmax=450,lb_db_ratio=3):
    '''
    calculates bed dimensions given v_total and nbv

    inputs
    v_total: volume to be treated [cm3]
    nbv: number of bed volumes to be treated [cm3(solute)/cm3(bed)]
    dbmax: maximum diameter of packed bed [cm], default = 450 cm
    lb_db_ratio: length to diameter ratio of packed bed [-], default = 3

    outputs
    db : diameter of packed bed [cm]
    lb : length of packed bed [cm]
    nbeds : number of beds [-]
    '''

    # define vmax
    vmax = (np.pi/4) * lb_db_ratio * dbmax**3

    # find volume of bed
    vb = v_total / nbv

    # calculate parameters for single bed systems
    nbeds = np.ones(vb.shape)

    db = ((4*vb)/(lb_db_ratio*np.pi))**(1/3)

    lb = lb_db_ratio*db

    # find indexes where vb > vmax (multiple beds needed)
    idx_overlim = vb > vmax

    # calculate parameters for multi-bed systems
    db[idx_overlim] = dbmax

    lb[idx_overlim] = dbmax * lb_db_ratio

    nbeds[idx_overlim] = vb[idx_overlim]/vmax

    return lb, db, nbeds
# lub_bed_dimensions()

def lmem_pore_diameter(Q, eps, n_a, v_s_bar, d_pore_mem_min_max=[None,None]):
    '''
    calculates the pore diameter of the adsorptive membrane based on saturation capacity
    (rearranged form of saturation_limit())

    inputs
    eps : membrane porosity [unitless*]
    n_a : Avogadro's number [1/mmol]
    v_s_bar : molar volume of solute [cm3/mmol]
    Q : saturation capacity of the adsorbent [mmol/cm3(particle)]
    d_pore_mem_min_max : lower and upper limits for pore diameter [cm], list floats, default [None(min),None(max)]

    outputs
    d_pore : diameter of pore [cm]
    '''
    d_pore = (4/Q) * (eps/(1-eps)) * (1/(n_a*np.pi))**(1/3) * (4/(3*v_s_bar))**(2/3) # [cm]

    # bound particle pore diameters if provided
    if d_pore_mem_min_max[0]:
        d_pore[d_pore<d_pore_mem_min_max[0]] = d_pore_mem_min_max[0]
        d_pore[d_pore>d_pore_mem_min_max[1]] = d_pore_mem_min_max[1]
    # End bounding particle pore diameter

    return d_pore
# lmem_pore_diameter()

def lmem_lem(N_BV,t_total,mu,dp,N_p,delP):
    '''
    calculates membrane thickness given process and design parameters

    inputs
    N_BV: vector of bed volumes [-], scalar or 1D numpy array
    t_total: regeneration time [s], scalar
    mu: fluid viscosity [Pa s], scalar
    dp : diameter of pore of membrane [cm]
    N_p: number of pores per unit area of membrane [1/cm^2], scalar
    delP: available pressure drop [Pa], scalar

    returns
    l_mem: membrane thickness [cm], scalar
    '''

    lmem = ((dp**4 * N_p * np.pi *  delP * t_total)/(128 * mu * N_BV))**(1/2) # [cm]

    return lmem
# END pore_size_limit()

def pb_bed_volume_limit(nbv,lub_bar,eps_b,eps_p,cin,v_s_bar,t_total,D_AB=1e-5,n_ps=1.46e9,n_a=6.022e20):
    '''
    return the residual of the bed volume limit equation for packed beds to use
    with a nonlinear equation solver
    nbv : number of bed volumes that can be treated [-*], float
    lub_bar : fractional length of unused bed, float
    eps_b : bed porosity [-*], float
    eps_p : particle porostiy [-*], float
    cin : inlet concentration of solute [mmol/cm3], float
    v_s_bar : molar volume of solute [cm3/mmol], float
    t_total : time for treatment [s], float
    D_AB: solute diffusivity [cm2/s], float, default 1e-5
    n_ps : number of pores per spherical adsorbent particle [-], default 1.46e9
    n_a : Avogadro's number [1/mmol], float, default 6.022e20

    returns
    res: residual of bed volume limit equation to solve with a nonlinear solver
    '''

    term1 = 4*eps_p*(1-eps_b)/cin
    term2 = ((3*n_ps*nbv)/(2*eps_p*lub_bar*D_AB*t_total))**(1/2)
    term3 = 1 / (n_a*np.pi)**(1/3)
    term4 = (4 / (3*v_s_bar))**(2/3)

    res = term1*term2*term3*term4 + eps_b - nbv

    return res
# pb_bed_volume_limit()

def crossover(eps, N_BV, r):
    '''
    Calculate the crossover point seen in Type 2 relative behavior between batch and semi-continuous processes.

    inputs
    eps: bed porosity, float scalar
    N_BV: number of bed volumes, float scalar
    r: removal ratio, float scalar
    '''

    xK = (eps - N_BV*(2-r))/(N_BV/r - eps);
    xQ = QB_target(xK, N_BV, r, epsilon=eps)

    return xK, xQ
# end crossover


## ------------------------------------ END model functions ------------------------------------ ##

## ---------------------------------------------------------------------------------------------- ##
##                                       workflow functions                                       ##
## ---------------------------------------------------------------------------------------------- ##

def calc_lithium_targets(cin_Li, v_s_bar_Li, delP_Li, l_mem, eps, m_mem_ratio_matrix,
t_total, mw_Li=6.941E-3, li_recovered=1000, batch_targets=False, removal_ratio=10,
n_points=2000, print_level=0, diagnostic_plots=False, calc_masses=False, path=None, casename=None):
    '''
    inputs
    cin_Li: concentration of lithium in feed [mmol/cm3], float scalar
    v_s_bar_Li: molar volume of lithium [cm3/mmol], float scalar
    delP_Li: pressure drop available for separation [Pa], float scalar
    l_mem: thickness of the membrane [m], float scalar
    eps: porosity of the memebrane [-], float scalar
    m_mem_ratio_matrix: quantity of membrane to be used for every unit of lithium recovered, integer numpy array, NOTE 2
    t_total: time to regeneration [s], float numpy array
    mw_Li: molecular mass of lithium [g/mmol], float scalar
    li_recovered: quantity of lithium to be recovered [kg], float scalar, NOTE 1
    batch_targets: flag to calculate batch material property targets, boolean, NOTE 5
    removal_ratio: ratio of inlet concentration to outlet concentration required for calculation of batch targets, boolean, NOTE 5
    n_points: number of points for sensitivity analysis contours, integer scalar
    print_level: verbosity of console outputs, integer scalar, NOTE 3
    diagnostic_plots: flag to print diagnostic plots, boolean
    calc_masses: calculate the masses of adsorbent needed for the given separation under consideration, bool, default False NOTE 6
    path: path to existing adsorbent data, string, default None, NOTE 6
    casename: unique name with which to save material requirement data, string, default None, NOTE 6

    returns
    K: vector of dimensionless binding affinity values used to calculate saturation capacity targets [-], float numpy array
    QC_Li: dimensionless saturation capacity targets  for a semicontinuous process corresponding to K [-], float numpy array
    qmax: upper bound of saturation capacity for the given system [mmol/g], float numpy array, NOTE 4
    K_qmax: vector of binding affinities corresponding to qmax [cm3/mmol], float numpy array
    QB_Li: dimensionless saturation capacity targets for a batch process corresponding to K [-], float numpy array, NOTE 5


    NOTES
    1. Amount of Li recovered: Based on the models, the results are insensitive to
    the amount of Li recovered. The user can convince themselves similarly by varying
    the value of li_recovered and noticing that it does not change the plots in any way.
    The insensitivity to amount of Li recovered is due to the value of N_BV staying
    constant which is a consequence of fixing the ratios of mass of Li recovered to mass
    of membrane used.

    2. mass of Li = total mass of Li recovered, not mass of Li recovered before regeneration

    3. 0 is lowest, 3 is highest for reasonable detail,
    4 will print unreasonably detailed information (think 1000 element floating point numpy arrays)

    4. system = membrane parameters + operating conditions

    5. QB_Li will be returned only if batch_targets == True

    6. Material requirement calculation was added during manuscript revisions.
    path and casename must be specified if calc_masses is set to True
    '''

    if print_level >= 3:
        print('Printing from calc_lithium_targets')
    # END print

    if print_level >= 2:
        print('System parameters for Li targets calculation:')

        disp_df = pd.DataFrame({
        'cin_Li [mmol/cm3]':cin_Li,
        'v_s_bar_Li [cm3/mmol]':v_s_bar_Li,
        'delP_Li [Pa]':delP_Li,
        'l_mem [m]':l_mem,
        'eps [-]':eps,
        't_total [s]':t_total})

        print(disp_df)
    # END print

    ## Get membrane parameters
    # Units: Np [1/m2], mu [Pa s], rho_mem [g/cm3], N_A [1/mmol], v_bar [cm3/mmol],
    # eps [-], A_mem [m2], sp_thr [m], cin_Pb [mmol/cm3], rho_mat [g/cm3]
    # Np, mu, rho_mem, N_A, __, __, __, __, __, rho_mat = get_membrane_model_params(get_rho_mat=True)
    Np = 1.6e14; # [1/m2] number of pores per m2 of membrane
    mu = 8.9e-4; # [Pa-s] viscosity of water at room temperature
    N_A = 6.022e20; # [1/mmol] Avogadro's number for 1 mmol of substance
    rho_mat = 1 # [g/cm3] (mat is the solid matrix of the polymer)

    # calculate sorbent density in CGS units
    rho_mem_updated = (1-eps)*rho_mat # [g/cm3]

    # convert units
    # rho_mem_SI = rho_mem * 1000 # [kg/m3]
    rho_mem_updated_SI = (1-eps)*rho_mat*1000 # [kg/m3]

    # vector of K values for sensitivity analysis calculations
    K = np.logspace(-2,8,n_points)
    nK = len(K)

    # calculate flow of brine required
    v_total = li_recovered * 1000 / (cin_Li*mw_Li) # [cm3]

    if print_level >= 3:
        print('li_recovered=',li_recovered,' kg')
        print('v_total=',v_total,' cm3')
    # END print

    # set ratio of mass of Li recovered to mass of membrane used
    # note that mass of Li = total mass of Li recovered, not mass of Li recovered before regeneration
    m_mem = m_mem_ratio_matrix*li_recovered # [kg]

    if print_level >= 3:
        print('m_mem [kg]=\n',m_mem)
    # END print

    ## volume of membrane
    # v_mem = m_mem  / rho_mem_SI # [m3]
    v_mem_updated = m_mem / rho_mem_updated_SI # [m3]

    if print_level >= 3:
        # print('v_mem [m3]=\n',v_mem)
        # print('rho_mem_SI [kg/m3] = \n',rho_mem_SI)
        print('v_mem_updated [m3]=\n',v_mem_updated)
        print('rho_mem_updated_SI [kg/m3] = \n',rho_mem_updated_SI)
    # END print

    # calculate nbv
    # nbv = v_total*1e-6/v_mem # [-]
    nbv_updated = v_total*1e-6/v_mem_updated # [-]

    if print_level >= 3:
        # print("nbv=\n",nbv)
        print('nbv_updated=\n',nbv_updated)
    # END print

    # preallocate storage
    nt = len(t_total)
    nnbv = len(nbv_updated)
    QC_Li = np.zeros((nK,nnbv))
    nbv_span_qmax = np.zeros((nt,nK))
    dp = np.zeros((nt,nK))
    qmax = np.zeros((nt,nK))
    K_qmax = np.zeros((nt,nK))

    # iterate over bed volumes and calculate dimensionless continuous material
    # property targets
    for j in range(nnbv):
        # QC_Li[:,j] = QC_target(K,nbv[j],epsilon=eps) # [-]
        QC_Li[:,j] = QC_target(K,nbv_updated[j],epsilon=eps) # [-]
    # END iteration over bed volumes

    # Set lower limit of nbv for qmax calculations based on epsilon value
    nbv_low_qmax = 1.7*eps # since nbv > eps is a requirement

    ## iterate over regeneration times and calculate Qmax
    for i in range(nt):

        ## calculate vector of N_BV for Qmax calculation
        # create lambda function to pass equation parameters
        nbv_limit = lambda root: bed_volume_limit(root, eps, cin_Li, t_total[i], Np, delP_Li, mu, l_mem, N_A, v_s_bar_Li)

        # solve nonlinear equaiton for upper limit of N_BV
        nbv_limit_soln = fsolve(nbv_limit,2000,full_output=True)

        if print_level >= 3:
            print('i=',i)
            print('nbv_limit_soln=',nbv_limit_soln)
        # END print

        # extract upper limit of N_BV
        nbv_high_qmax = 0.99999*nbv_limit_soln[0][0]

        # generate vector of N_BV
        nbv_span_qmax[i,:] = np.linspace(nbv_low_qmax,nbv_high_qmax,n_points)

        # calculate dp
        dp[i,:] = pore_size_limit(nbv_span_qmax[i,:],t_total[i],mu,l_mem,Np,delP_Li) # [m]

        if np.amin(dp[i,:] <= 5e-9):
            print('WARNING!!!\nPore size smaller than 5 nm for t_total[i] = {0} [s]'.format(t_total[i]))

            fig_dp, ax_dp = plt.subplots()
            plt.plot(nbv_span_qmax[i,:], dp[i,:])
            plt.xlabel('N_BV [-]')
            plt.ylabel('dp [m]')
            ax_dp.set_yticks([1e-9,3e-9,5e-9,7e-9])
            plt.ylim(1e-9,7e-9)
            plt.grid(True)
            plt.plot()
        # END print

        if diagnostic_plots:
            print('nothing to show')
        # END diagnostic plots

        # calculate Qmax
        qmax[i,:] = saturation_limit(eps,dp[i,:],N_A,v_s_bar_Li) # [mmol/cm3]

        # calculate K
        K_qmax[i,:] = calc_K(nbv_span_qmax[i,:],qmax[i,:],eps,cin_Li) # [cm3/mmol]

        # save arrays to get input values for dry run calculations
        np.savetxt(casename+'qmax'+str(i)+'.csv',qmax[i,:],delimiter=',')
        np.savetxt(casename+'K_qmax'+str(i)+'.csv',K_qmax[i,:],delimiter=',')
        np.savetxt(casename+'nbv_span_qmax'+str(i)+'.csv',nbv_span_qmax[i,:],delimiter=',')
        np.savetxt(casename+'nbv_span_by_t_qmax_'+str(i)+'.csv',nbv_span_qmax[i,:]/t_total[i],delimiter=',')
        np.savetxt(casename+'dp'+str(i)+'.csv',dp[i,:],delimiter=',')

        if print_level >= 4:
            print('i=',i)
            print('K_qmax=\n',K_qmax)
        # END print
    # END iterate over regeneration times and calculate Qmax

    # convert qmax to mmol/g
    qmax = qmax / rho_mat # [mmol/g]

    # calculate material requirements
    if calc_masses:
        ## Existing lithium adsorbents
        # read existing Li adsorbent data
        li_adsorbent_data = pd.read_csv(path)
        # separate into meaningful variables
        li_Ks = li_adsorbent_data['K [cm3/mmol_Li]'].to_numpy() # cm3/mmol
        li_Qs = li_adsorbent_data['Q [mmol_Li/cm3_ads]'].to_numpy() # mmol/cm3
        li_legend = li_adsorbent_data['Legend'].tolist() # -

        m_ads = calc_semicont_masses(v_total, rho_mem_updated, eps, li_Qs, li_Ks, cin_Li)

        pd.DataFrame({'Legend':li_legend,
        'K[cm3/mmol]':li_Ks,
        'Q[mmol/cm3]':li_Qs,
        'm_ads[g]':m_ads}).to_csv(path_or_buf='semicont_material_requirements_'+casename+'.csv')
        print('Saved material mass requirements to semicont_material_requirements_'+casename+'.csv')
        # np.savetxt('semicont_material_requirements_'+casename+'.csv',
        # [li_legend,li_Ks,li_Qs,m_ads],
        # header="Legend,K[cm3/mmol],Q[mmol/cm3],m_ads[g]",
        # delimiter=',')
    # END calc_masses

    if batch_targets:

        # preallocate storage
        QB_Li = np.zeros((nK,nnbv))

        # iterate over bed volumes and calculate dimensionless batch material
        # property targets
        for j in range(nnbv):
            QB_Li[:,j] = QB_target(K, nbv_updated[j], removal_ratio, epsilon=eps) # [dimensionless]
        # END iteration over bed volumes

        return K, QC_Li, K_qmax, qmax, QB_Li
    else:
        return K, QC_Li, K_qmax, qmax
    # END return after calculating batch targets if required

# END calc_lithium_targets()

def calc_qmax(casename,v_total,v_s_bar,c_in,
              n_el,t_total,
              Dp_min=0.03,Dp_max=0.15,
              print_level=0,design_bed=False,detailed_plots=False,lub_bar=0.01):
    '''
    calculates Qmax of the packed bed given design and operating criteria

    inputs
    v_total :  volume to be treated [cm3]
    v_s_bar : molar volume of solute [cm3/mmol]
    c_in: inlet concentration of solute [mmol/cm3]
    '''

    # define constants for this problem
    n_ps = 1.46e9 # [-], number of pores per spherical particle
    mu = 8.9e-4 # [Pa-s], viscosity
    eps_part = 0.57 # [-*], porosity of particle
    eps_bed = 0.4 # [-*], bed void fraction
    n_a = 6.022e20 # [1/mmol] Avogadro's number
    D_AB = 1e-5 # [cm2/s] solute diffusivity

    nbv_min = bed_volumes(D_AB,t_total,Dp_max,lub_bar=lub_bar)
    nbv_max = bed_volumes(D_AB,t_total,Dp_min,lub_bar=lub_bar)

    if nbv_min < eps_bed:
        # maintain nbv > eps constraint
        nbv_min = 1.1*eps_bed

    if print_level >= 3:
        print('N_BV_min=',nbv_min)
        print('N_BV_max=',nbv_max)
    # End print

    n_bv = np.linspace(nbv_min,nbv_max,n_el)

    ## find particle diameter corresponding to n_bv
    d_particle = particle_diameter(D_AB, t_total, n_bv, lub_bar=lub_bar)

    if print_level >= 3:
        print('particle diameter [cm]=',d_particle)
    # end print

    ## find pore diameter
    d_pore = particle_pore_diameter(eps_part, n_ps, d_particle) # [cm]

    # mask array indexes where pore diameter becomes smaller than 5 nm
    d_pore_lower_than_5 = d_pore < 5e-7
    # save to file for plotting
    np.savetxt('./Pb_packed_bed_qmax_underlimit_pore'+casename+'.csv',d_pore_lower_than_5,delimiter=',',fmt='%d')

    if print_level >= 3:
        print('particle pore diameter [cm]=',d_pore)
    # End print

    ## find Qmax
    Qmax = saturation_capacity(eps_part, d_pore, n_a, v_s_bar)

    if print_level >= 3:
        print('saturation capacity [mmol/cm3(particle)]=',Qmax)
    # End print


    ## find K(Qmax)
    Kmax = binding_affinity(n_bv, eps_bed, Qmax, c_in)

    if print_level >= 3:
        print('Kmax [cm3(solvent)/mmol] = ',Kmax)
    # End

    # save to file for plotting
    np.savetxt('./Pb_packed_bed_qmax'+casename+'.csv',np.array([Qmax,Kmax,n_bv]),delimiter=',')

    ## design bed
    if design_bed:

        # find bed dimensions
        db, lb, n_bed= bed_dimensions(v_total,n_bv,lb_db_ratio=3)

        if print_level >= 3:
            print('db [cm] = ',db)
            print('lb [cm] = ',lb)
            print('n_bed [-] = ',n_bed)
        # end print

        # find pressure drop in the tower
        delP = pressure_drop(lb, mu, eps_bed, n_bv, d_particle, t_total)

        # mask array indexes where pressure drop exceeds 30 psi
        delP_more_than_30 = delP/6894.76 > 30

        # set flag to check if any pressure drops exceed 30 psi
        if delP_more_than_30.max() > 0:
            delP_exceeds_30 = True
        else:
            delP_exceeds_30 = False
        # end

        # save to file for plotting
        np.savetxt('./Pb_packed_bed_qmax_overlimit_pressure'+casename+'.csv',delP_more_than_30,delimiter=',',fmt='%d')

        if print_level >= 3:
            print('delP [psi]= ',delP/6894.76)
        # End print

        # find Peclet number in the tower
        n_pe = peclet_number(d_particle,v_total,n_bed,D_AB,db,t_total)

        if print_level >= 3:
            print('Pe = ',n_pe)
        # End print

        if print_level >= 2:
            print('when pore size is under the limit of 5nm =')
            print('saturation capacity [mmol/cm3(particle)]=',Qmax[d_pore_lower_than_5][0])
            print('Kmax [cm3(solvent)/mmol] = ',Kmax[d_pore_lower_than_5][0])
            print('N_BV = ',n_bv[d_pore_lower_than_5][0])
            print('db [cm] = ',db[d_pore_lower_than_5][0])
            print('lb [cm] = ',lb[d_pore_lower_than_5][0])
            print('n_bed [-] = ',n_bed[d_pore_lower_than_5][0])
            print('delP [psi]= ',delP[d_pore_lower_than_5][0]/6894.76)
            print('Pe = ',n_pe[d_pore_lower_than_5][0])

            if delP_exceeds_30:
                print('when pressure drop in bed is more than 30 psi =')
                print('particle size [mm]=',d_particle[delP_more_than_30][0]*10)
                print('pore size [nm]=',d_pore[delP_more_than_30][0]*1e7)
                print('saturation capacity [mmol/cm3(particle)]=',Qmax[delP_more_than_30][0])
                print('Kmax [cm3(solvent)/mmol] = ',Kmax[delP_more_than_30][0])
                print('N_BV = ',n_bv[delP_more_than_30][0])
                print('db [cm] = ',db[delP_more_than_30][0])
                print('lb [cm] = ',lb[delP_more_than_30][0])
                print('n_bed [-] = ',n_bed[delP_more_than_30][0])
                print('delP [psi]= ',delP[delP_more_than_30][0]/6894.76)
                print('Pe = ',n_pe[delP_more_than_30][0])
            # end
        # End print

        if detailed_plots:

            # Q_bar vs K_bar
            plt.figure(figsize=(10,10))
            plt.loglog(Kmax*c_in, Qmax/c_in)
            plt.loglog(Kmax[d_pore_lower_than_5][0]*c_in, Qmax[d_pore_lower_than_5][0]/c_in,
                       'kx', markersize=12, label='pore under 5nm')

            if delP_exceeds_30:
                plt.loglog(Kmax[delP_more_than_30][0]*c_in, Qmax[delP_more_than_30][0]/c_in,
                           'ro', markersize=12, label='delP over 30psi')
            # end

            plt.xlabel('K_bar',fontsize=24)
            plt.ylabel('Q_bar',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()

             # K_bar vs n_bv
            plt.figure(figsize=(10,10))
            plt.semilogy(n_bv, Kmax*c_in)
            plt.semilogy(n_bv[d_pore_lower_than_5][0], Kmax[d_pore_lower_than_5][0]*c_in,
                        'kx', markersize=12, label='pore under 5nm')
            if delP_exceeds_30:
                plt.semilogy(n_bv[delP_more_than_30][0], Kmax[delP_more_than_30][0]*c_in,
                           'ro', markersize=12, label='delP over 30psi')
            # end
            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('K_bar',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()

            # n_bed vs n_bv
            plt.figure(figsize=(10,10))
            plt.plot(n_bv, n_bed,'o')
            plt.plot(n_bv[d_pore_lower_than_5][0], n_bed[d_pore_lower_than_5][0],
                    'kx', markersize=12, label='pore under 5nm')
            if delP_exceeds_30:
                plt.plot(n_bv[delP_more_than_30][0], n_bed[delP_more_than_30][0],
                       'ro', markersize=12, label='delP over 30psi')
            # end
            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('n_bed',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()


            # pore diameter vs n_bv
            plt.figure(figsize=(10,10))
            plt.plot(n_bv, d_pore*1e7)
            plt.plot(n_bv[d_pore_lower_than_5][0], d_pore[d_pore_lower_than_5][0]*1e7,
                    'kx', markersize=12, label='pore under 5nm')
            if delP_exceeds_30:
                plt.plot(n_bv[delP_more_than_30][0], d_pore[delP_more_than_30][0]*1e7,
                       'ro', markersize=12, label='delP over 30psi')
            # end
            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('Pore dia. (nm)',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()

            # particle diameter vs n_bv
            plt.figure(figsize=(10,10))
            plt.plot(n_bv, d_particle*10)
            plt.plot(n_bv[d_pore_lower_than_5][0], d_particle[d_pore_lower_than_5][0]*10,
                    'kx', markersize=12, label='pore under 5nm')

            if delP_exceeds_30:
                plt.plot(n_bv[delP_more_than_30][0], d_particle[delP_more_than_30][0]*10,
                       'ro', markersize=12, label='delP over 30psi')
            # end

            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('Particle dia. (mm)',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()

            # peclet number vs n_bv
            plt.figure(figsize=(10,10))
            plt.plot(n_bv, n_pe)
            plt.plot(n_bv[d_pore_lower_than_5][0], n_pe[d_pore_lower_than_5][0],
                    'kx', markersize=12, label='pore under 5nm')

            if delP_exceeds_30:
                plt.plot(n_bv[delP_more_than_30][0],  n_pe[delP_more_than_30][0],
                       'ro', markersize=12, label='delP over 30psi')
            # end

            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('Peclet number',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()

            # delta P vs n_bv
            plt.figure(figsize=(10,10))
            plt.plot(n_bv, delP/6894.76)
            plt.plot(n_bv[d_pore_lower_than_5][0], delP[d_pore_lower_than_5][0]/6894.76,
                    'kx', markersize=12, label='pore under 5nm')
            if delP_exceeds_30:
                plt.plot(n_bv[delP_more_than_30][0],  delP[delP_more_than_30][0]/6894.76,
                       'ro', markersize=12, label='delP over 30psi')
            # End

            plt.xlabel('n_bv',fontsize=24)
            plt.ylabel('delP (psi)',fontsize=24)
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.grid(True,which='both')
            plt.show()
# end calc_qmax()

def plot_lithium_targets(K, QC_Li, K_qmax, qmax, m_mem_ratio_matrix, cin_Li, eps, t_total,
path, title, plot_batch=False, QB_Li=None, ls = ['--', ':','-.', (0, (1, 10))],
w_in=3., h_in=4.0, dpi_fig=1200, print_level=0, dimensionless=False,pbed_qmax_path=None,
pore_underlimit_path=None,delP_overlimit_path=None,pb_labels=None,combined_plot=False,
mem_Qmax_labels=None):
    '''
    inputs
    K: vector of dimensionless binding affinity values used to calculate saturation capacity targets [-], float numpy array
    QC_Li: dimensionless saturation capacity targets for a semicontinuous process corresponding to K [-], float numpy array
    K_qmax: vector of binding affinities corresponding to qmax [cm3/mmol], float numpy array
    qmax: upper bound of saturation capacity for the given system [mmol/g], float numpy array, NOTE 1
    m_mem_ratio_matrix: quantity of membrane to be used for every unit of lithium recovered, integer numpy array, NOTE 3
    cin_Li: concentration of lithium in feed [mmol/cm3], float scalar
    eps: Porosity of the membrane [-], float scalar
    t_total: time to regeneration [s], float numpy array
    path: path at which existing lithium adsorbent data is saved, string
    title: title for the plot, also used as filename to save plot, string
    plot_batch: flag to plot batch targets. boolean
    QB_Li: dimensionless saturation capacity targets for a batch process corresponding to K [-], float numpy array
    ls: vector of linestyles for Qmax lines, list
    w_in: figure width [inches], float scalar
    h_in: figure height [inches], float scalar
    dpi_fig: pixel density [dpi], float scalar
    print_level: verbosity of console outputs, integer scalar, NOTE 2
    dimensionless: Flag to toggle dimensionless plots on and off

    returns
    fig, ax: matplotlib figure and axis objects with plots of lithium recovery targets with dimensions


    NOTES
    1. system = membrane parameters + operating conditions

    2. 0 is lowest, 3 is highest for reasonable detail,
    4 will print unreasonably detailed information (think 1000 element floating point numpy arrays)

    3. mass of Li = total mass of Li recovered, not mass of Li recovered before regeneration
    '''

    if print_level >= 3:
        print('Printing from plot_lithium_targets()')
    # END print

    if print_level >= 4:

        fig_d0, ax_d0 = plt.subplots()
        plt.loglog(K,QC_Li)
        plt.xlabel('K [-]')
        plt.ylabel('QC_Li [-]')
        plt.grid(True)
        plt.savefig(title+'-argument_ref_1.png',dpi=300,bbox_inches='tight')
        plt.show()

        fig_d1, ax_d1 = plt.subplots()
        plt.loglog(K_qmax,qmax)
        plt.xlabel('K_qmax [cm3/mmol]')
        plt.ylabel('qmax [mmol/g]')
        plt.grid(True)
        plt.savefig(title+'-argument_ref_2.png',dpi=300,bbox_inches='tight')
        plt.show()

        print('m_mem_ratio_matrix = \n',m_mem_ratio_matrix)

        disp_df = pd.DataFrame({'cin_Li [mmol/cm3]':cin_Li,
        'eps [-]':eps,
        't_total [s]':t_total,
        'path':path,
        'title':title,
        'plot_batch':plot_batch})

        print(disp_df)
    # END print

    ## Get membrane parameters
    # Units: Np [1/m2], mu [Pa s], rho_mem [g/cm3], N_A [1/mmol], v_bar [cm3/mmol],
    # eps [-], A_mem [m2], sp_thr [m], cin_Pb [mmol/cm3], rho_mat [g/cm3]
    # __, __, __, __, __, __, __, __, __, rho_mat = get_membrane_model_params(get_rho_mat=True)
    rho_mat = 1 # [g/cm3] (mat is the solid matrix of the polymer)

    # find updated membrane density based on porosity and matrix density
    rho_mem = (1-eps)*rho_mat # [g/cm3]

    ## Existing lithium adsorbents
    # read existing Li adsorbent data
    li_adsorbent_data = pd.read_csv(path)
    # separate into meaningful variables
    li_Ks = li_adsorbent_data['K [cm3/mmol_Li]'].to_numpy() # cm3/mmol
    li_Qs = li_adsorbent_data['Q [mmol_Li/cm3_ads]'].to_numpy() # mmol/cm3
    # convert Q to mmol/g
    li_Qs = li_Qs/rho_mat # [mmol/g]

    # convert to dimensionless variables if required
    if dimensionless:
        li_Ks = li_Ks*cin_Li # [-]
        li_Qs = li_Qs*rho_mat/cin_Li # [-]

        if print_level >= 3:
            print("li_Ks:",li_Ks)
            print("li_Qs:",li_Qs)
        # END
    # END conversion to dimensionless variables

    li_legend = li_adsorbent_data['Legend'].tolist() # -

    # create matplotlib objects
    fig, ax = plt.subplots(figsize=(w_in,h_in),dpi=dpi_fig)

    # loop over bed volumes and plot contours
    for j in range(len(m_mem_ratio_matrix)):

        # make legend entries pretty
        mpt_ratio = '1:{0:d}'.format(m_mem_ratio_matrix[j])

        # convert from dimensionless to with units if required
        if not dimensionless:
            Kplot = K/cin_Li # [cm3/mmol]

            QCplot = QC_Li[:,j]*cin_Li/rho_mat # [mmol/g]
        else:
            Kplot = K # [-]

            QCplot = QC_Li[:,j] # [-]
        # END convert from dimensionless to with units if required

        # plot property targets
        line = plt.loglog(Kplot, QCplot, label=mpt_ratio)

        if plot_batch:

            # convert from dimensionless to with units if required
            if not dimensionless:
                QBplot = QB_Li[:,j]*cin_Li/rho_mat # [mmol/g]
            else:
                QBplot = QB_Li[:,j] # [-]
            # END convert from dimensionless to with units if required

            plt.loglog(Kplot,QBplot,'--',color=line[-1].get_color())
        # END plot batch targets if required

    # END loop over bed volumes and plot contours

    if combined_plot:

        # plot Packed Bed contours
        if pbed_qmax_path:

            # pb_linestyles = ['r--','b--','m--']
            # fig_pb, ax_pb = plt.subplots(figsize=(w_in,h_in),dpi=dpi_fig)
            # plt.xscale('log')
            # plt.yscale('log')

            # read contour data
            K_x = np.loadtxt(pbed_qmax_path+'K_x.csv',delimiter=',')
            Q_y = np.loadtxt(pbed_qmax_path+'Q_y.csv',delimiter=',')
            lub_bar = np.loadtxt(pbed_qmax_path+'lub_bar.csv',delimiter=',')

            # convert to dimensionless coordinates if required
            if dimensionless:
                K_x = K_x*cin_Li
                Q_y = Q_y/cin_Li
            # end dimensionless conversion

            # plot
            plt.contour(K_x,Q_y,lub_bar,levels=[1e-2,1e-1],colors=['b','b'],linestyles=['--',':'])
            plt.plot([],[],'b--',label=r'$\overline{l_{ub}}$=1%')
            plt.plot([],[],'b:',label=r'$\overline{l_{ub}}$=10%')
            # plt.plot
        else:
            # plot membrane contours
            # loop over membrane thicknesses and plot Qmax lines
            for i in range(len(mem_Qmax_labels)):

                # convert time into nice legend entries
                # if mem_Qmax_labels[i] == '1 mm':
                #     label_fmt = r'$l_b$=1 mm'
                # elif mem_Qmax_labels[i] == '0.3 mm':
                #     label_fmt = r'$l_b$=0.3 mm'
                # # End nice legend
                label_fmt = r'$l_b$='+mem_Qmax_labels[i]

                # convert from dimensionless to with units if required
                if not dimensionless:
                    qmax_plot = qmax[i,:] # [mmol/g]

                    K_qmax_plot = K_qmax[i,:] # [cm3/mmol]
                else:
                    qmax_plot = qmax[i,:]*rho_mat/cin_Li # [-]

                    K_qmax_plot = K_qmax[i,:]*cin_Li # [-]

                    # save to file to verify with dry runs
                    np.savetxt('qmax_plot'+str(i),qmax_plot,delimiter=',')
                    np.savetxt('K_qmax_plot'+str(i),K_qmax_plot,delimiter=',')
                # END convert from dimensionless to with units if required

                # plot
                plt.loglog(K_qmax_plot,qmax_plot,'k',linestyle=ls[i],label=label_fmt)
                # plt.loglog(K_qmax[i,:],qmax[i,:],'k',linestyle=ls[i],label='Regen every {0} {1}'.format(time_val,time_unit))
            # loop over membrane thicknesses and plot Qmax lines
    else:
        # loop over regeneration times and plot Qmax lines
        for i in range(len(t_total)):

            # convert time into nice legend entries
            if (np.round(t_total[i]/3600) < 1):
                time_val = int(np.round(t_total[i]/60))
                time_unit = 'min(s)'
            elif (np.round(t_total[i]/3600) < 24):
                time_val = int(np.round(t_total[i]/3600))
                time_unit = 'hour(s)'
            elif ((np.round(t_total[i]/3600/24) <= 7)):
                time_val = int(np.round(t_total[i]/3600/24))
                time_unit = 'day(s)'
            elif ((np.round(t_total[i]/3600/24/30) < 12)):
                time_val = int(np.round(t_total[i]/3600/24/30))
                time_unit = 'month(s)'
            elif((np.round(t_total[i]/3600/24/30) >= 12)):
                time_val = int(np.round(t_total[i]/3600/24/30/12))
                time_unit = 'year(s)'
            # End nice legend

            # more legend formatting
            label_1 = r'$t_{bt}$ = '
            label_2 = '{0} {1}'.format(time_val,time_unit)
            label = label_1+label_2

            # convert from dimensionless to with units if required
            if not dimensionless:
                qmax_plot = qmax[i,:] # [mmol/g]

                K_qmax_plot = K_qmax[i,:] # [cm3/mmol]
            else:
                qmax_plot = qmax[i,:]*rho_mat/cin_Li # [-]

                K_qmax_plot = K_qmax[i,:]*cin_Li # [-]

                # save to file to verify with dry runs
                np.savetxt('qmax_plot'+str(i),qmax_plot,delimiter=',')
                np.savetxt('K_qmax_plot'+str(i),K_qmax_plot,delimiter=',')
            # END convert from dimensionless to with units if required

            # plot
            plt.loglog(K_qmax_plot,qmax_plot,'k',linestyle=ls[i],label=label)
            # plt.loglog(K_qmax[i,:],qmax[i,:],'k',linestyle=ls[i],label='Regen every {0} {1}'.format(time_val,time_unit))
        # loop over regeneration times and plot Qmax lines
    # end plotting qmax lines

    # plot Li adsorbent data
    for k in range(len(li_legend)):
        # letters of the alphabet markers
        # plt.loglog(li_Ks[k],li_Qs[k],'k',marker='${0}$'.format(li_legend[k]),markersize=8)

        # dots as markers
        plt.loglog(li_Ks[k],li_Qs[k],'kx',markersize=3)
    # END plot Li adsorbent data

    # # plot Packed Bed contours
    # if pbed_qmax_path:
    #
    #     # pb_linestyles = ['r--','b--','m--']
    #     # fig_pb, ax_pb = plt.subplots(figsize=(w_in,h_in),dpi=dpi_fig)
    #     # plt.xscale('log')
    #     # plt.yscale('log')
    #
    #     # read contour data
    #     K_x = np.loadtxt(pbed_qmax_path+'K_x.csv',delimiter=',')
    #     Q_y = np.loadtxt(pbed_qmax_path+'Q_y.csv',delimiter=',')
    #     lub_bar = np.loadtxt(pbed_qmax_path+'lub_bar.csv',delimiter=',')
    #
    #     # convert to dimensionless coordinates if required
    #     if dimensionless:
    #         K_x = K_x*cin_Li
    #         Q_y = Q_y/cin_Li
    #     # end dimensionless conversion
    #
    #     # plot
    #     plt.contour(K_x,Q_y,lub_bar,levels=[1e-2,1e-1],colors=['b','b'],linestyles=['--',':'])
    #     plt.plot([],[],'b--',label=r'$\overline{l_{ub}}$=1%')
    #     plt.plot([],[],'b:',label=r'$\overline{l_{ub}}$=10%')
    #     # plt.plot
    #
    # # plot Packed Bed contours

    # make plots pretty

    # choose axis labels based on type of plot
    if not dimensionless:
        plot_xlabel = r"$\mathbf{K}$ (cm$\mathbf{^3}$ mmol$\mathbf{^{-1}}$)"
        plot_ylabel = r"$\mathbf{Q}$ (mmol / g$\mathbf{^{-1}}$)"
    else:
        plot_xlabel = r"$\mathbf{\overline{K}}$"
        plot_ylabel = r"$\mathbf{\overline{Q}}$"
    # choose axis labels based on type of plot

    plt.xlabel(plot_xlabel,fontsize=12,fontweight='bold')
    plt.ylabel(plot_ylabel,fontsize=12,fontweight='bold')

    # square aspect ratio
    set_aspect_ratio_log(ax, 1.0)

    # ticklabels and axis limits
    # plt.xlim(1e-2,1e7)
    # plt.ylim(1e-2,1e3)
    # ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
    # ax.set_xticks([1e-2,1e0,1e2,1e4,1e6,1e8])

    # legend
    # ncol=3
    ncol = 2

    ## row major legend (print across columns)
    handles, labels = ax.get_legend_handles_labels()
    # handles = flip(handles, ncol)
    # labels = flip(labels, ncol)

    # formatting for 3 column legened
    # plt.legend(handles, labels, loc='upper center',borderaxespad=0.,ncol=ncol,
    # fontsize=10, bbox_to_anchor=(0.35,-0.23),
    # handletextpad=0.5,columnspacing=0.5,framealpha=1.0,title='Mass of Li recovered:Mass of membrane')

    # formatting for 2 column legend
    plt.legend(handles, labels, loc='upper center',borderaxespad=0.,ncol=ncol,
    fontsize=10, bbox_to_anchor=(0.37,-0.23),
    handletextpad=0.5,columnspacing=0.5,framealpha=1.0)

    plt.grid(True,which='major')
    # plt.grid(True,which='minor',linewidth=0.5)
    plt.minorticks_off()

    plt.savefig(title+'.png',bbox_inches='tight')

    return fig, ax
# END plot_lithium_targets()

def make_log_contourplot(x,y,z,xlabel,ylabel,zlabel,outlines=None,plot_existing=False,
                         K_exist=None,Q_exist=None,save=None):

    fig,ax = plt.subplots(figsize=(10,10))
    plt.xscale('log')
    plt.yscale('log')

    cs = plt.contourf(x,y,z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)

    if outlines:
        contour = ax.contour(x,y,z,levels=outlines,colors='black')
        # plt.clabel(contour,fontsize=18,colors='k',fmt='%.e')
    # end

    cbar = fig.colorbar(cs)
    cbar.set_label(label = zlabel, size=24)
    cbar.ax.tick_params(labelsize=24)

    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlabel(xlabel,fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.ylim(1e-2,5e1)

    if plot_existing:
        # Placeholder for sorbents in legend
        ph = list(map(chr, range(65,90)))

        for i in range(len(K_exist)):
            # letters of the alphabet as markers
            plt.loglog(K_exist[i],Q_exist[i],"k",marker="${0}$".format(ph[i]),markersize=14,)
#             plt.loglog(K_exist[i],Q_exist[i],"k",marker="x",markersize=14,)
        # end loop over existing sorbents
    # end plotting existing sorbents


    if save:
        plt.savefig(save+'.jpg',dpi=300,bbox_inches='tight')
    # end save

    plt.show()

    return fig, ax
# end make_log_contourplot()

def qmax_heatmaps(casename,eps_bed,eps_part,eps_mem,c_in,v_s_bar,
                v_total,t_total,delP_mem,
                D_AB=1e-5,mu=8.9e-4,n_ps=1.46e9,Np_mem=1.6e10,
                n_a=6.022e20,K_lim=[-2,8],Q_lim=[-3,3],n_el=100,plot_existing=False,
                K_exist=None,Q_exist=None,
                d_pore_min_max=[None,None],
                d_pore_mem_min_max=[None,None],
                plot_contours=True):
    '''
    plots heatmaps of qmax and intermediates for packed bed and membrane adsorbents
    by sampling K and Q

    inputs
    casename : name of the case being run, used as filename prefix, string
    eps_bed : bed porosity [-*], float
    eps_part : particle porostiy [-*], float
    eps_mem : membrane_porosity [-*], float
    c_in : inlet concentration of solute [mmol/cm3], float
    v_s_bar : molar volume of solute [cm3/mmol], float
    v_total : volume to be treated [cm3], float
    t_total : time for treatment [s], float
    delP_mem : pressure drop available for membrane separations [Pa], float
    D_AB: solute diffusivity [cm2/s], float, default 1e-5
    mu : viscosity of water at room temperature [Pa s], float, default 8.9e-4
    n_ps : number of pores per spherical adsorbent particle [-], default 1.46e9
    Np_mem : number of pores per unit area of membrane [1/cm2], float, default 1.6e10
    n_a : Avogadro's number [1/mmol], float, default 6.022e20
    K_lim: 10^upper and 10^lower limits for binding affinity (K) grid, list of ints, default [-2,8] (based on Pb study)
    Q_lim: 10^upper and 10^lower limits for saturation capacity (Q) grid, list of ints, default [-3,3] (based on Pb study)
    n_el : number of elements for grid vector generation, int, default 100
    plot_existing : flag to plot existing sorbent data on heatmaps, bool, default False
    K_exist : K values of existing sorbents, np.array of floats
    Q_exist : Q values of existing sorbents, np.array of floats
    d_pore_min_max : lower and upper limits for pore diameter of particle [cm], list floats, default [None(min),None(max)]
    d_pore_mem_min_max : lower and upper limits for pore diameter of membrane [cm], list floats, default [None(min),None(max)]
    plot_contours: flag to plot contour maps of calculated values, bool, default True
    '''

    # sample K, Q on log space
    K_arr = np.logspace(K_lim[0],K_lim[1],n_el)
    Q_arr = np.logspace(Q_lim[0],Q_lim[1],n_el)

    # create meshgrid
    K_x, Q_y = np.meshgrid(K_arr, Q_arr)
    # save to file
    np.savetxt(casename+'K_x.csv',K_x,delimiter=',')
    print('Saved K values to '+casename+'K_x.csv')
    np.savetxt(casename+'Q_y.csv',Q_y,delimiter=',')
    print('Saved Q values to '+casename+'Q_y.csv')

    # calculate and plot N_BV given K, Q
    nbv = nbv_semicont(K_x,Q_y,eps_bed,c_in) # -
    if plot_contours:
        fig_nbv, ax_nbv = make_log_contourplot(K_x,Q_y,nbv,'K [1/M]','Q [cm3/mmol]','N_BV',
                                           outlines= [1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6],
                                           plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                           save=casename+'nbv_map')


   # calculate dp given Q
    dp = lub_pore_diameter(Q_y,eps_part,n_a,v_s_bar,d_pore_min_max=d_pore_min_max) # cm
    if plot_contours:
        fig_dp, ax_dp = make_log_contourplot(K_x,Q_y,dp*1e7,'K [1/M]','Q [cm3/mmol]','dp_PB [nm]',
                                         outlines=[1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'dpore_PB_map')

    # calculate Dp given dp
    Dp = lub_particle_diameter(eps_part,n_ps,dp) # cm
    if plot_contours:
        fig_Dp, ax_Dp = make_log_contourplot(K_x,Q_y,Dp*10,'K [1/M]','Q [cm3/mmol]','Dp [mm]',
                                         outlines=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'Dpart_map')

    # find lub_bar
    lub_bar = lub_lub_bar(Dp,nbv,D_AB,t_total) # -
    if plot_contours:
        fig_lub, ax_lub = make_log_contourplot(K_x,Q_y,lub_bar,
                                         r'K [M$^{-1}$]',r'Q [mmol cm$^{-3}$]',r'$\overline{l_{ub}}$',
                                         outlines=[1e-4,1e-3,1e-2,1e-1,1e0,1e1],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'lub_map')

    # save to file
    np.savetxt(casename+'lub_bar.csv',lub_bar,delimiter=',')
    print('Saved lub_bar values to '+casename+'lub_bar.csv')

    # find bed dimensions
    lb, db, nbeds = lub_bed_dimensions(v_total,nbv) # cm, cm, -

    # find pressure drop
    delP = pressure_drop(lb,mu,eps_bed,nbv,Dp,t_total) # Pa
    if plot_contours:
        fig_delP, ax_delP = make_log_contourplot(K_x,Q_y,delP/6894.76,'K [1/M]','Q [cm3/mmol]','delP [psi]',
                                         outlines=[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'delP_PB_map')

    # find pore diameter of membrane
    dp_mem = lmem_pore_diameter(Q_y, eps_mem, n_a, v_s_bar, d_pore_mem_min_max=d_pore_mem_min_max) # cm
    if plot_contours:
        fig_dp_mem, ax_dp_mem = make_log_contourplot(K_x,Q_y,dp_mem*1e7,'K [1/M]','Q [cm3/mmol]','dp_mem [nm]',
                                         outlines=[1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'dpore_MEM_map')

    # find thickness of membrane
    lmem = lmem_lem(nbv,t_total,mu,dp_mem,Np_mem,delP_mem) # cm
    if plot_contours:
        fig_lmem, ax_lmem = make_log_contourplot(K_x,Q_y,lmem*10,'K [1/M]','Q [cm3/mmol]','lmem [mm]',
                                         outlines=[1e-2,1e-1,1e0,1e1],
                                         plot_existing=plot_existing,K_exist=K_exist,Q_exist=Q_exist,
                                         save=casename+'l_MEM_map')

# end qmax_heatmaps()


def calc_lead_targets(cin_Pb, v_s_bar_Pb, delP_Pb, l_mem, eps, m_mem,
v_total, t_total, batch_targets=False, removal_ratio=None,
n_points=2000, print_level=0, calc_masses=False,
path=None, casename=None):
    '''
    inputs
    cin_Pb: concentration of lead in feed [mmol/cm3], float scalar
    v_s_bar_Pb: molar volume of lead [cm3/mmol], float scalar
    delP_Li: pressure drop available for separation [Pa], float scalar
    l_mem: thickness of the membrane [m], float scalar
    eps: porosity of the memebrane [-], float scalar
    m_mem: quantity of membrane to be used for separation, integer numpy array
    v_total: volume of water to be treated [cm3], float scalar
    t_total: time to regeneration [s], float numpy array
    batch_targets: flag to calculate batch material property targets, boolean, NOTE 1
    removal_ratio: ratio of inlet concentration to outlet concentration required for calculation of batch targets, boolean, NOTE 1
    n_points: number of points for sensitivity analysis contours, integer scalar
    print_level: verbosity of console outputs, integer scalar, NOTE 2
    calc_masses: calculate the masses of adsorbent needed for the given separation under consideration, bool, default False NOTE 3
    path: path to existing adsorbent data, string, default None, NOTE 3
    casename: unique name with which to save material requirement data, string, default None, NOTE 3

    returns
    K: vector of dimensionless binding affinity values used to calculate saturation capacity targets [-], float numpy array
    QC_Pb: dimensionless saturation capacity targets  for a semicontinuous process corresponding to K [-], float numpy array
    qmax: upper bound of saturation capacity for the given system [mmol/g], float numpy array, NOTE 4
    K_qmax: vector of binding affinities corresponding to qmax [cm3/mmol], float numpy array
    QB_Pb: dimensionless saturation capacity targets for a batch process corresponding to K [-], float numpy array, NOTE 1


    NOTES
    1. QB_Li will be returned only if batch_targets == True

    2. 0 is lowest, 3 is highest for reasonable detail,
    4 will print unreasonably detailed information (think 1000 element floating point numpy arrays)

    3. removal_ratio, path and casename must be specified if calc_masses is set to True
    '''

    if print_level >= 3:
        print('Printing from calc_lead_targets')
    # END print

    if print_level >= 2:
        print('System parameters for Pb targets calculation:')

        disp_df = pd.DataFrame({
        'cin_Pb [mmol/cm3]':cin_Pb,
        'v_s_bar_Pb [cm3/mmol]':v_s_bar_Pb,
        'delP_Pb [Pa]':delP_Pb,
        'l_mem [m]':l_mem,
        'eps [-]':eps,
        'v_total [m3]':v_total,
        't_total [s]':t_total})

        print(disp_df)
    # END print

    ## Get membrane parameters
    # Units: Np [1/m2], mu [Pa s], rho_mem [g/cm3], N_A [1/mmol], v_bar [cm3/mmol],
    # eps [-], A_mem [m2], sp_thr [m], cin_Pb [mmol/cm3], rho_mat [g/cm3]
    # Np, mu, rho_mem, N_A, __, __, __, __, __, rho_mat = get_membrane_model_params(get_rho_mat=True)

    # Define constants for the problem
    Np = 1.6e14; # [1/m2] number of pores per m2 of membrane
    mu = 8.9e-4; # [Pa-s] viscosity of water at room temperature
    N_A = 6.022e20; # [1/mmol] Avogadro's number for 1 mmol of substance
    rho_mat = 1 # [g/cm3] (mat is the solid matrix of the polymer)

    # calculate sorbent density in CGS units
    rho_mem_updated = (1-eps)*rho_mat # [g/cm3]

    # convert units
    # rho_mem_SI = rho_mem * 1000 # [kg/m3]
    rho_mem_updated_SI = (1-eps)*rho_mat*1000 # [kg/m3]

    # vector of K values for sensitivity analysis calculations
    K = np.logspace(-4,8,n_points)
    nK = len(K)

    if print_level >= 3:
        print('m_mem [kg]=\n',m_mem)
    # END print

    ## volume of membrane
    # v_mem = m_mem  / rho_mem_SI # [m3]
    v_mem_updated = m_mem / rho_mem_updated_SI # [m3]

    if print_level >= 3:
        # print('v_mem [m3]=\n',v_mem)
        # print('rho_mem_SI [kg/m3] = \n',rho_mem_SI)
        print('v_mem_updated [m3]=\n',v_mem_updated)
        print('rho_mem_updated_SI [kg/m3] = \n',rho_mem_updated_SI)
    # END print

    # calculate nbv
    # nbv = v_total*1e-6/v_mem # [-]
    nbv_updated = v_total/v_mem_updated # [-]

    if print_level >= 3:
        # print("nbv=\n",nbv)
        print('nbv_updated=\n',nbv_updated)
    # END print

    # preallocate storage
    nt = len(t_total)
    nnbv = len(nbv_updated)
    QC_Pb = np.zeros((nK,nnbv))
    nbv_span_qmax = np.zeros((nt,nK))
    dp = np.zeros((nt,nK))
    qmax = np.zeros((nt,nK))
    K_qmax = np.zeros((nt,nK))

    # iterate over bed volumes and calculate dimensionless continuous material
    # property targets
    for j in range(nnbv):
        QC_Pb[:,j] = QC_target(K,nbv_updated[j],epsilon=eps) # [-]
    # END iteration over bed volumes

    # Set lower limit of nbv for qmax calculations based on epsilon value
    nbv_low_qmax = 1.7*eps # since nbv > eps is a requirement

    ## iterate over regeneration times and calculate Qmax
    for i in range(nt):

        ## calculate vector of N_BV for Qmax calculation
        # create lambda function to pass equation parameters
        nbv_limit = lambda root: bed_volume_limit(root, eps, cin_Pb, t_total[i], Np, delP_Pb, mu, l_mem, N_A, v_s_bar_Pb)

        # solve nonlinear equaiton for upper limit of N_BV
        nbv_limit_soln = fsolve(nbv_limit,2000,full_output=True)

        if print_level >= 3:
            print('i=',i)
            print('nbv_limit_soln=',nbv_limit_soln)
        # END print

        # extract upper limit of N_BV
        nbv_high_qmax = 0.99999*nbv_limit_soln[0][0]

        # generate vector of N_BV
        nbv_span_qmax[i,:] = np.linspace(nbv_low_qmax,nbv_high_qmax,n_points)

        # calculate dp
        dp[i,:] = pore_size_limit(nbv_span_qmax[i,:],t_total[i],mu,l_mem,Np,delP_Pb) # [m]

        if np.amin(dp[i,:] <= 5e-9):
            print('WARNING!!!\nPore size smaller than 5 nm for t_total[i] = {0} [s]'.format(t_total[i]))

            fig_dp, ax_dp = plt.subplots()
            plt.plot(nbv_span_qmax[i,:], dp[i,:])
            plt.xlabel('N_BV [-]')
            plt.ylabel('dp [m]')
            ax_dp.set_yticks([1e-9,3e-9,5e-9,7e-9])
            plt.ylim(1e-9,7e-9)
            plt.grid(True)
            plt.plot()
        # END print

        # calculate Qmax
        qmax[i,:] = saturation_limit(eps,dp[i,:],N_A,v_s_bar_Pb) # [mmol/cm3]

        # calculate K
        K_qmax[i,:] = calc_K(nbv_span_qmax[i,:],qmax[i,:],eps,cin_Pb) # [cm3/mmol]

        if print_level >= 4:
            print('i=',i)
            print('K_qmax=\n',K_qmax)
        # END print
    # END iterate over regeneration times and calculate Qmax

    # convert qmax to mmol/g
    qmax = qmax / rho_mat # [mmol/g]

    # calculate material requirements
    if calc_masses:
        ## Existing lead adsorbents
        # read list of sorbents
        # read sorbent property data from a .csv file
        sorbent_dat = pd.read_csv(path,sep=',').values.tolist()

        # extract data to meaningful variable names
        sorbents = [];
        sorb_Qs = [];
        sorb_Ks = [];
        sorb_refs = [];

        for i in range(len(sorbent_dat)):
            sorbents.append(sorbent_dat[i][0]);
            sorb_Qs.append(sorbent_dat[i][1]); # mmol/g
            sorb_Ks.append(sorbent_dat[i][2]); # cm3/mmol
            sorb_refs.append(sorbent_dat[i][3]);
        # END read sorbent data

        # convert lists to numpy arrays
        Pb_sorb_Qs = np.array(sorb_Qs,dtype=float)
        Pb_sorb_Ks = np.array(sorb_Ks,dtype=float)

        # convert Q to volume basis
        Pb_sorb_Qs = Pb_sorb_Qs*rho_mat # mmol/cm3

        m_ads = calc_semicont_masses(v_total*1e6, rho_mem_updated, eps, Pb_sorb_Qs, Pb_sorb_Ks, cin_Pb)

        pd.DataFrame({'K[cm3/mmol]':Pb_sorb_Ks,
        'Q[mmol/cm3]':Pb_sorb_Qs,
        'm_ads[g]':m_ads}).to_csv(path_or_buf='semicont_material_requirements_'+casename+'.csv')
        print('Saved semicontinuous material requirements to '+'semicont_material_requirements_'+casename+'.csv')

        m_batch = calc_batch_masses(v_total*1e6,  rho_mem_updated, eps, Pb_sorb_Qs, Pb_sorb_Ks, cin_Pb, cin_Pb/removal_ratio)
        pd.DataFrame({'K[cm3/mmol]':Pb_sorb_Ks,
        'Q[mmol/cm3]':Pb_sorb_Qs,
        'm_batch[g]':m_batch}).to_csv(path_or_buf='batch_material_requirements_'+casename+'.csv')
        print('Saved batch material requirements to '+'batch_material_requirements_'+casename+'.csv')
    # END calc_masses

    if batch_targets:

        # preallocate storage
        QB_Pb = np.zeros((nK,nnbv))

        # iterate over bed volumes and calculate dimensionless batch material
        # property targets
        for j in range(nnbv):
            QB_Pb[:,j] = QB_target(K, nbv_updated[j], removal_ratio, epsilon=eps) # [dimensionless]
        # END iteration over bed volumes

        return K, QC_Pb, K_qmax, qmax, QB_Pb
    else:
        return K, QC_Pb, K_qmax, qmax
    # END return after calculating batch targets if required

# END calc_lithium_targets()

def plot_lead_targets(K, QC_Pb, K_qmax, qmax, m_mem_matrix, cin_Pb, eps, v_total,
t_total, path, title, rho_mat = 1, plot_batch=False, QB_Pb=None, ls = ['--', ':','-.', (0, (1, 10))],
w_in=3., h_in=4.0, dpi_fig=1200, print_level=0, dimensionless=False, pbed_qmax_path=None,
pb_labels=None,combined_plot=False, mem_Qmax_labels=None):
    '''
    inputs
    K: vector of dimensionless binding affinity values used to calculate saturation capacity targets [-], float numpy array
    QC_Pb: dimensionless saturation capacity targets for a semicontinuous process corresponding to K [-], float numpy array
    K_qmax: vector of binding affinities corresponding to qmax [cm3/mmol], float numpy array
    qmax: upper bound of saturation capacity for the given system [mmol/g], float numpy array
    m_mem_matrix: quantity of membrane to be used for separation, integer numpy array
    cin_Pb: concentration of lithium in feed [mmol/cm3], float scalar
    eps: Porosity of the membrane [-], float scalar
    v_total: volume of water to be treated [m3], float scalar
    t_total: time to regeneration [s], float numpy array
    path: path at which existing lead adsorbent data is saved, string
    title: title for the plot, also used as filename to save plot, string
    plot_batch: flag to plot batch targets. boolean
    QB_Li: dimensionless saturation capacity targets for a batch process corresponding to K [-], float numpy array
    ls: vector of linestyles for Qmax lines, list
    w_in: figure width [inches], float scalar
    h_in: figure height [inches], float scalar
    dpi_fig: pixel density [dpi], float scalar
    print_level: verbosity of console outputs, integer scalar
    dimensionless: Flag to toggle dimensionless plots on and off

    returns
    fig, ax: matplotlib figure and axis objects with plots of lithium recovery targets with dimensions


    NOTES

    '''

    if print_level >= 3:
        print('Printing from plot_lead_targets()')
    # END print

    if print_level >= 4:

        fig_d0, ax_d0 = plt.subplots()
        plt.loglog(K,QC_Pb)
        plt.xlabel('K [-]')
        plt.ylabel('QC_Pb [-]')
        plt.grid(True)
        plt.savefig(title+'-argument_ref_1.png',dpi=300,bbox_inches='tight')
        plt.show()

        fig_d1, ax_d1 = plt.subplots()
        plt.loglog(K_qmax,qmax)
        plt.xlabel('K_qmax [cm3/mmol]')
        plt.ylabel('qmax [mmol/g]')
        plt.grid(True)
        plt.savefig(title+'-argument_ref_2.png',dpi=300,bbox_inches='tight')
        plt.show()

        print('m_mem_matrix = \n',m_mem_matrix)

        disp_df = pd.DataFrame({'cin_Pb [mmol/cm3]':cin_Pb,
        'eps [-]':eps,
        't_total [s]':t_total,
        'v_total [m3]':v_total,
        'path':path,
        'title':title,
        'plot_batch':plot_batch})

        print(disp_df)
    # END print

    # find updated membrane density based on porosity and matrix density
    rho_mem = (1-eps)*rho_mat # [g/cm3]

    ## Existing lead adsorbents
    # read list of sorbents
    # read sorbent property data from a .csv file
    sorbent_dat = pd.read_csv(path,sep=',').values.tolist()

    # extract data to meaningful variable names
    sorbents = [];
    sorb_Qs = [];
    sorb_Ks = [];
    sorb_refs = [];

    for i in range(len(sorbent_dat)):
        sorbents.append(sorbent_dat[i][0]);
        sorb_Qs.append(sorbent_dat[i][1]); # mmol/g
        sorb_Ks.append(sorbent_dat[i][2]); # cm3/mmol
        sorb_refs.append(sorbent_dat[i][3]);
    # END read sorbent data

    # convert lists to numpy arrays
    Pb_sorb_Qs = np.array(sorb_Qs,dtype=float)
    Pb_sorb_Ks = np.array(sorb_Ks,dtype=float)

    # convert to dimensionless if required
    if dimensionless:
        Pb_sorb_Qs = Pb_sorb_Qs*rho_mat/cin_Pb
        Pb_sorb_Ks = Pb_sorb_Ks*cin_Pb
    # dimensionless conversion

    # create matplotlib objects
    fig, ax = plt.subplots(figsize=(w_in,h_in),dpi=dpi_fig)

    # loop over bed volumes and plot contours
    for j in range(len(m_mem_matrix)):

        # make legend entries pretty
        label = '{0:d} kg'.format(m_mem_matrix[j])

        # convert from dimensionless to with units if required
        if not dimensionless:
            Kplot = K/cin_Pb # [cm3/mmol]

            QCplot = QC_Pb[:,j]*cin_Pb/rho_mat # [mmol/g]
        else:
            Kplot = K # [-]

            QCplot = QC_Pb[:,j] # [-]
        # END convert from dimensionless to with units if required

        # plot property targets
        line = plt.loglog(Kplot, QCplot, label=label)

        if plot_batch:

            # convert from dimensionless to with units if required
            if not dimensionless:
                QBplot = QB_Pb[:,j]*cin_Pb/rho_mat # [mmol/g]
            else:
                QBplot = QB_Pb[:,j] # [-]
            # END convert from dimensionless to with units if required

            plt.loglog(Kplot,QBplot,'--',color=line[-1].get_color())
        # END plot batch targets if required

    # END loop over bed volumes and plot contours

    if pbed_qmax_path:
        # plot Packed Bed contours
        # read contour data
        K_x = np.loadtxt(pbed_qmax_path+'K_x.csv',delimiter=',')
        Q_y = np.loadtxt(pbed_qmax_path+'Q_y.csv',delimiter=',')
        lub_bar = np.loadtxt(pbed_qmax_path+'lub_bar.csv',delimiter=',')

        # convert to dimensionless coordinates if required
        if dimensionless:
            K_x = K_x*cin_Li
            Q_y = Q_y/cin_Li
        # end dimensionless conversion

        # plot
        plt.contour(K_x,Q_y,lub_bar,levels=[1e-2,1e-1],colors=['b','b'],linestyles=['--',':'])
        plt.plot([],[],'b--',label=r'$\overline{l_{ub}}$=1%')
        plt.plot([],[],'b:',label=r'$\overline{l_{ub}}$=10%')
        # plt.plot
    else:
        # plot membrane contours
        # loop over membrane thicknesses and plot Qmax lines
        for i in range(len(mem_Qmax_labels)):

            # # convert time into nice legend entries
            # if mem_Qmax_labels[i] == '1 mm':
            #     label_fmt = r'$l_b$=1 mm'
            # elif mem_Qmax_labels[i] == '0.3 mm':
            #     label_fmt = r'$l_b$=0.3 mm'
            # elif
            # # End nice legend
            label_fmt = r'$l_b$='+mem_Qmax_labels[i]

            # convert from dimensionless to with units if required
            if not dimensionless:
                qmax_plot = qmax[i,:] # [mmol/g]

                K_qmax_plot = K_qmax[i,:] # [cm3/mmol]
            else:
                qmax_plot = qmax[i,:]*rho_mat/cin_Li # [-]

                K_qmax_plot = K_qmax[i,:]*cin_Li # [-]
            # END convert from dimensionless to with units if required

            # plot
            plt.loglog(K_qmax_plot,qmax_plot,'k',linestyle=ls[i],label=label_fmt)
            # plt.loglog(K_qmax[i,:],qmax[i,:],'k',linestyle=ls[i],label='Regen every {0} {1}'.format(time_val,time_unit))
        # loop over membrane thicknesses and plot Qmax lines
    # plot membrane contours

    # plot lead adsorbent data
    for k in range(len(Pb_sorb_Ks)):
        # letters of the alphabet markers
        # plt.loglog(li_Ks[k],li_Qs[k],'k',marker='${0}$'.format(li_legend[k]),markersize=8)

        # dots as markers
        plt.loglog(Pb_sorb_Ks[k],Pb_sorb_Qs[k],'kx',markersize=3)
    # END plot Li adsorbent data

    # make plots pretty

    # choose axis labels based on type of plot
    if not dimensionless:
        plot_xlabel = r"$\mathbf{K}$ (M$\mathbf{^{-1}}$)"
        plot_ylabel = r"$\mathbf{Q}$ (mmol g$\mathbf{^{-1}}$)"
    else:
        plot_xlabel = r"$\mathbf{\overline{K}}$"
        plot_ylabel = r"$\mathbf{\overline{Q}}$"
    # choose axis labels based on type of plot

    plt.xlabel(plot_xlabel,fontsize=12,fontweight='bold')
    plt.ylabel(plot_ylabel,fontsize=12,fontweight='bold')

    # square aspect ratio
    set_aspect_ratio_log(ax, 1.0)

    # legend
    # ncol=3
    ncol = 2

    ## row major legend (print across columns)
    handles, labels = ax.get_legend_handles_labels()

    # formatting for 2 column legend
    plt.legend(handles, labels, loc='upper center',borderaxespad=0.,ncol=ncol,
    fontsize=10, bbox_to_anchor=(0.37,-0.23),
    handletextpad=0.5,columnspacing=0.5,framealpha=1.0)

    plt.grid(True,which='major')
    # plt.grid(True,which='minor',linewidth=0.5)
    plt.minorticks_off()

    plt.savefig(title+'.png',bbox_inches='tight')

    return fig, ax
# END plot_lead_targets()


def nondim_sens(rec,eps,nbv,K):
    '''
    sensitivity analysis for dimensionless materials property targets

    inputs
    rec: removal ratio, float scalar
    eps: bed porosity, float scalar
    nbv: number of bed volumes to calculate targets over, float numpy array
    K: binding affinity vector for saturation capacity (Q) targets calculation, float numpy array

    outputs
    QC: dimensionless saturation capacity targets  for a semicontinuous process corresponding to K [-], float numpy array
    QB: dimensionless saturation capacity targets for a batch process corresponding to K, float numpy array
    xKQ: [K,Q] value at which crossover occurs, float numpy array
    '''

    # preallocate storage
    nK = len(K)
    nnbv = len(nbv)
    QC = np.zeros((nK,nnbv))
    QB = np.zeros((nK,nnbv))
    xKQ = np.zeros((2,nnbv))

    for j in range(nnbv):
        # semicontinuous targets
        QC[:,j] = QC_target(K,nbv[j],epsilon=eps) # [-]

        # batch targets
        QB[:,j] = QB_target(K, nbv[j], rec, epsilon=eps) # [dimensionless]

        xKQ[0,j], xKQ[1,j] = crossover(eps, nbv[j], rec)
    # END iteration over bed volumes

    return QC, QB, xKQ
# end nondim_sens()

def plot_dimensionless_targets(K, QC, QB, xKQ, nbv, rec, eps, w_in=3., h_in=4.0, save=None):
    '''
    make plots for the dimensionless material property target contours in
    the supporting information

    inputs
    K: vector of dimensionless binding affinity values used to calculate saturation capacity targets [-], float numpy array
    QC: dimensionless saturation capacity targets for a semicontinuous process corresponding to K [-], float numpy array
    QB: dimensionless saturation capacity targets for a batch process corresponding to K [-], float numpy array
    xKQ: [K,Q] value at which crossover occurs, float numpy array
    nbv: number of bed volumes to calculate targets over, float numpy array
    rec: removal ratio, float scalar
    eps: bed porosity, float scalar
    w_in: figure width [inches], float scalar, default 3
    h_in: figure height [inches], float scalar, default 4
    save: filename to save plot, string, default None (don't save)

    outputs
    dimensionless material property target contours, saved to .png file if
    filename is provided in the kwargs
    '''

    # create matplotlib objects
    fig, ax = plt.subplots(figsize=(w_in,h_in))

    # line colors
    colors = ['#1f77b4', '#ff7f0e', '#8c564b', '#e377c2',\
              '#069740', '#1e0027', '#2ca02c', '#d62728']

    # loop over bed volumes and plot contours
    for j in range(len(nbv)):
        # legend label
        mpt_label = r'$N_{BV}=$'+str(nbv[j])

        # plot semicontinuous property targets
        line = plt.loglog(K, QC[:,j], label=mpt_label, color=colors[j])

        # plot batch property targets
        plt.loglog(K,QB[:,j],'--',color=line[-1].get_color())

        # crossover point
        plt.loglog(xKQ[0,j],xKQ[1,j],marker="o",\
                   markeredgecolor=line[-1].get_color(),markerfacecolor="none",markersize=5);

    # end loop over bed volumes and plot contours

    # axes labels and title
    plot_xlabel = r"$\mathbf{\overline{K}}$"
    plot_ylabel = r"$\mathbf{\overline{Q}}$"
    plt.xlabel(plot_xlabel,fontsize=12,fontweight='bold')
    plt.ylabel(plot_ylabel,fontsize=12,fontweight='bold')
    plt_title = str('r='+str(rec)+r' $\mathbf{\epsilon}$='+str(eps))
    plt.title(plt_title,fontsize=16,fontweight='bold')

    # grid and axis ticks
    x_lim_high = 5*1e4
    x_lim_low = 5*1e-5
    ax.set_xticks([1e-4,1e-2,1e0,1e2,1e4]) # x-ticks
    ax.set_xlim(x_lim_low,x_lim_high) # x-lim
    ax.tick_params(labelsize=12) # ticksize
    plt.grid(True,which='major')
    plt.minorticks_off()

    ## legend
    # Ghost points for legend
    plt.plot([],[],'k--',label='Batch')
    plt.plot([],[],'k-',label='Semi-\ncontinuous')

    # if crossover is plotted, include it in the legend
    if (x_lim_low < np.amax(xKQ[0,:]) and\
        np.amax(xKQ[0,:]) < x_lim_high):
        plt.plot([],[],marker="o",markeredgecolor='k',markerfacecolor="none",\
                 markersize=5,linestyle='none',label='Predicted\nCrossover')
    # END: if crossover is plotted, include it in the legend

    # print legend row-first
    # set number of columns in the legend
    ncol = 3
    # get handles and labels of the legend to flip
    handles, labels = ax.get_legend_handles_labels()
    # plot legend
    plt.legend(flip(handles,ncol),flip(labels,ncol),\
               loc='upper center',ncol=ncol,fontsize=10,borderaxespad=0.,\
               handletextpad=0.5,columnspacing=0.5,\
               bbox_to_anchor = (0.35,-0.25),
               labelspacing=0.2,handlelength=1)

    # square plots
    set_aspect_ratio_log(ax, 1.0)

    if save:
        plt.savefig(save+'.png',bbox_inches='tight',dpi=300)
    # end save to file

    # display on screen
    plt.show()

    return fig, ax
## ------------------------------------ END workflow functions ------------------------------------ ##
