#=============================================================================================
# Tools for producing gromacs .mdp files for free energy calculations.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-01-22
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.utilities import Units
from mmtools.moltools.ligandtools import *
import commands
import random

random.seed()

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================



#=============================================================================================
# SUBROUTINES
#=============================================================================================

def writeFile(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()

    return

def mdpWriter(parameters):
    """Produce a gromacs .mdp file for free energy calculations.

    TODO
      Use Units.py for things like dt, duration.
      Autodetect comm_mode based on vacuum or solvent.
      Set seed randomly to unique seed.
    """

    # parameter block
    parameters = dict()    

    # set default parameters
    parameters['title'] = ""
    parameters['minimization_steps'] = 1000
    parameters['dt'] = 0.002 * Units.fs
    parameters['duration'] = 2.0 * Units.ns
    parameters['comm_mode'] = "Linear"
    parameters['seed'] = random.randint(0,2**32-1)
    parameters['temperature'] = 298.0 * Units.K
    parameters['pressure'] = 1.0 * Units.atm 

    # compuare various properties in gromacs units
    parameters['nsteps'] = int(parameters['duration'] / parameters['dt'])     
    parameters['dt_in_ps'] = parameters['dt'] / Units.ps
    parameters['duration_in_ns'] = parameters['duration'] / Units.ns
    parameters['temperature_in_K'] = parameters['temperature'] / Units.K
    parameters['pressure_in_bar'] = parameters['pressure'] / Units.bar

    # Construct header.
    header_block = """
; VARIOUS PREPROCESSING OPTIONS = 
title                    = %(title)s
cpp                      = /usr/bin/cpp
define                   =
""" 

    # Construct minimizer block.
    minimizer_block = """
; RUN CONTROL PARAMETERS = 
integrator               = steep
nsteps                   = %(minimization_steps)s

; ENERGY MINIMIZATION OPTIONS = 
; Force tolerance and initial step-size = 
emtol                    = 10
emstep                   = 1.0e-4
; Max number of iterations in relax_shells = 
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints = 
fcstep                   = 0
; Frequency of steepest descents steps when doing CG = 
nstcgsteep               = 1000
""" 
    
    # dynamics control block
    nsteps = int(round(duration_in_ns * 1000.0 / dt)) # compute number of steps
    dynamics_block = """
; RUN CONTROL PARAMETERS = 
integrator               = vv
; start time and timestep in ps = 
tinit                    = 0
dt                       = %(dt_in_fs)f ; fs
nsteps                   = %(nsteps)d ; %(duration_in_ns).3f ns
; mode for center of mass motion removal = 
comm-mode                = %(comm_mode)s
; number of steps for center of mass motion removal = 
nstcomm                  = 500
; group(s) for center of mass motion removal = 
comm-grps                = 
""" 

    # output control block
    output_block = """
; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 500
nstenergy                = 500
; Output frequency and precision for xtc file = 
nstxtcout                = 50
xtc-precision            = 10000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
xtc-grps                 = System
; Selection of energy groups = 
energygrps               = System
""" 
    
    # thermal and pressure control block
    bath_block = """
; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = %(temperature_in_K)f ; K
gen_seed                 = %(seed)d

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = AndersenII
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 5.0
ref_t                    = %(temperature_in_K)f ; K

; Pressure coupling      = 
Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 1.67 ; ps
compressibility          = 4.5e-5 ; 1/bar
ref_p                    = %(pressure_in_bar)f ; bar

; LANGEVIN DYNAMICS OPTIONS = 
; Temperature, friction coefficient (amu/ps) and random seed = 
bd-temp                  = %(temperature_in_K)f
bd-fric                  = 0
""" 

    # nonbonded parameters
    nonbonded_block_solvated = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 0.9
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 0.75
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = switch
; cut-off lengths        = 
rvdw-switch              = 0.7
rvdw                     = 0.75
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = AllEnerPres
; Spacing for the PME/PPPM FFT grid = 
fourierspacing           = 0.10
; FFT grid size, when a value is 0 fourierspacing will be used = 
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters = 
pme_order                = 6
ewald_rtol               = 1e-06
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = no
"""    

    nonbonded_block_vacuum = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10 ; can we set this to 0 for vacuum?
; ns algorithm (simple or grid) = 
ns_type                  = simple
; Periodic boundary conditions: xyz or no = 
pbc                      = no
; nblist cut-off         = 
rlist                    = 25.0 ; to emulate no cutoff
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = cut-off
rcoulomb-switch          = 0
rcoulomb                 = 23.0 ; to emulate no cutoff
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = cut-off
; cut-off lengths        = 
rvdw-switch              = 0
rvdw                     = 23.0 ; to emulate no cutoff
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = no
"""    
    # free energy block
    free_energy_block = """
; OPTIONS FOR EXPANDED ENSEMBLE SIMULATIONS
; Free energy control stuff = 
free-energy              = mutate ; annihilate electrostatics and Lennard-Jones
nstfep                   = 50 ; 0.1 ps between weight updates (must be integer multiple of nstlist)
nstdgdl                  = 50 ; 0.1 ps between writing energies

; weight update scheme
lambda-mc                = gibbs-wang-landau ; Wang-Landau with waste recycling
mc-wldelta               = 1.0 ; initial delta factor for Wang-Landau
mc-wlscale               = 0.5 ; amount by which delta is adjusted for Wang-Landau
mc-nratio                = 0.8 ; flatness criterion

; state transition probability
move-mc                  = metropolized-gibbs ; Metropolized Gibbs for fastest mixing of states

; starting and stopping
mc-nstart                = 0 ; number of updates to perform per state for driving through each state
mc-nequil                = 1000 ; number of updates before freezing weights

init-lambda              = 1 ; initial state

; schedule for switching off lambdas
; first, restraints are turned on as charges are switched off
; next, vdw and torsions are switched off
fep-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
coul-lambda              = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
restraint-lambda         = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
vdw-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
bonded-lambda            = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

sc-alpha                 = 0.5 ; soft-core factor
"""
    
    # construct .mdp parameter files
    unconstrained_minimization = (header_block + minimize_block + output_block + nonbonded_block) % parameters
    writeFile('minimize-unconstrained.mdp', unconstrained_minimization)
    constrained_minimization = (header_block + minimize_block + output_block + nonbonded_block + constraint_block) % parameters
    writeFile('minimize-constrained.mdp', constrained_minimization)    
    production = (header_block + dynamics_block + output_block + nonbonded_block + constraint_block + free_energy_block) % parameters        
    writeFile('production.mdp', production)

    return


