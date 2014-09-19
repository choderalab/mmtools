#=============================================================================================
# Example illustrating use of 'moltools' in computing absolute hydration free energies of small molecules.
#
# PROTOCOL
#
# * Construct molecule from IUPAC name (protonation and tautomeric states are heuristically guessed) [OpenEye OEMol]
# * Generate multiple likely conformations to use for replicates [OpenEye Omega]
# * Parameterize the molecule with GAFF [AmberTools Antechamber]
# * Set up the solvated system [AmberTools LEaP]
# * Convert the system to gromacs topology and coordinates
# * Set up vacuum and solvated runs with gromacs_dg
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-01-23
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.moltools.ligandtools import *
from mmtools.gromacstools.MdpFile import *    
import commands
import os    
from numpy import *

#=============================================================================================
# CHANGELOG
#============================================================================================
# 2008-04-19 JDC
# * Both vacuum and solvated simulations now have counterions added for charged solutes.
#=============================================================================================

#=============================================================================================
# KNOWN BUGS
#============================================================================================
# * A vacuum simulation should not be set up for decoupling simulations.
# * Multiple solute conformations should be used to see different replicates and the solute should start in the decoupled state.
#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

# gromacs executables
mdrun = 'mdrun_d'
grompp = 'grompp_d'
editconf = 'editconf_d'

nreplicates = 3 # number of replicates of each calculation to set up

clearance = 5.0 # clearance around solute for box construction, in Angstroms

#=============================================================================================
# GROMACS PARAMETERS
#=============================================================================================

# Header block
mdpblock = dict()

mdpblock['header']  = """
; VARIOUS PREPROCESSING OPTIONS = 
title                    = title
cpp                      = /usr/bin/cpp
define                   =
""" 

# Construct minimizer block.
mdpblock['minimization'] = """
; RUN CONTROL PARAMETERS = 
integrator               = steep
nsteps                   = 1000

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
mdpblock['dynamics'] = """
; RUN CONTROL PARAMETERS = 
integrator               = vv
; start time and timestep in ps = 
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000 ; 5 ns
; mode for center of mass motion removal = 
comm-mode                = Linear
; number of steps for center of mass motion removal = 
nstcomm                  = 500
; group(s) for center of mass motion removal = 
comm-grps                = system
""" 

# output control block
mdpblock['output'] = """
; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 500 ; (must be integral multiple of nstdgdl for post-analysis)
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 500 ; 
nstenergy                = 50 ; (must be same as nstdgdl for analysis scripts)
; Output frequency and precision for xtc file = 
nstxtcout                = 50 ; 0.1 ps
xtc-precision            = 1000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
xtc-grps                 = solute
; Selection of energy groups = 
energygrps               = System
""" 

# constraints block
mdpblock['constraints'] = """
; OPTIONS FOR BONDS     = 
constraints              = hbonds ; constrain bonds to hydrogen
; Type of constraint algorithm = 
constraint-algorithm     = shake
; Do not constrain the start configuration = 
unconstrained-start      = no
; Use successive overrelaxation to reduce the number of shake iterations = 
Shake-SOR                = no
; Relative tolerance of shake = 
shake-tol                = 1e-12
; Highest order in the expansion of the constraint coupling matrix = 
lincs-order              = 4
; Lincs will write a warning to the stderr if in one step a bond = 
; rotates over more degrees than = 
lincs-warnangle          = 180
; Convert harmonic bonds to morse potentials = 
morse                    = no
"""
    
# thermal and pressure control block
mdpblock['thermostat'] = """
; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300.0 ; K
gen_seed                 = 12345

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = AndersenII
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 1.0
ref_t                    = 300.0 ; K
"""

# thermal and pressure control block
mdpblock['thermostat-equilibration'] = """
; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300.0 ; K
gen_seed                 = 12345

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = AndersenII
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 1.0
ref_t                    = 300.0 ; K
"""

mdpblock['barostat'] = """
; Pressure coupling      = 
Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 1.67 ; ps
compressibility          = 4.5e-5 ; 1/bar
ref_p                    = 1.0 ; bar
""" 

# dmobley protocol
mdpblock['nonbonded-solvent'] = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.0 ; MRS has 0.9
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 0.9 ; MRS has 0.75
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = switch
; cut-off lengths        = 
rvdw-switch              = 0.8 ; MRS has 0.7
rvdw                     = 0.9 ; MRS has 0.75
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

# MRS protocol
mdpblock['nonbonded-solvent-mrs'] = """
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

# dmobley protocol
mdpblock['nonbonded-vacuum'] = """
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10 ; can we set this to 0 for vacuum?
; ns algorithm (simple or grid) = 
ns_type                  = simple
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 25.0 ; to emulate no cutoff
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = cut-off
rcoulomb-switch          = 0
rcoulomb                 = 25.0 ; to emulate no cutoff
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon-r                = 1
; Method for doing Van der Waals = 
vdw-type                 = cut-off
; cut-off lengths        = 
rvdw                     = 25.0
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = no
"""    


# JDC protocol
mdpblock['nonbonded-vacuum-jdc'] = """
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
mdpblock['free-energy'] = """
; OPTIONS FOR EXPANDED ENSEMBLE SIMULATIONS
; Free energy control stuff = 
free-energy              = mutate ; annihilate electrostatics and Lennard-Jones
nstfep                   = 50 ; 0.1 ps between weight updates (must be integer multiple of nstlist)
nstdgdl                  = 50 ; 0.1 ps between writing energies (must be same as nstdgdl for analysis scripts)

; weight update scheme
lambda-mc                = gibbs-wang-landau ; Wang-Landau with waste recycling
mc-wldelta               = 1.0 ; initial delta factor for Wang-Landau (in kT)
mc-wlscale               = 0.5 ; scalar by which delta is scaled for Wang-Landau
mc-nratio                = 0.8 ; flatness criterion -- histograms are reset after all states are sampled within mc-nratio factor of the mean

; state transition probability
move-mc                  = metropolized-gibbs ; Metropolized Gibbs for fastest mixing of states

; starting and stopping
mc-nstart                = 0 ; number of updates to perform per state for driving through each state
mc-nequil                = 2000000 ; number of steps before freezing weights (4 ns)

init-lambda              = 1 ; initial state

; schedule for switching off lambdas
; first, restraints are turned on as charges are switched off
; next, vdw and torsions are switched off
;fep-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
;coul-lambda              = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
;restraint-lambda         = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
;vdw-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
;bonded-lambda            = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

; more Coulomb
fep-lambda               = 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
coul-lambda              = 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
restraint-lambda         = 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
vdw-lambda               = 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
bonded-lambda            = 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

; works okay with annihilation of charged benzamidine
;fep-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
;coul-lambda              = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
;restraint-lambda         = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
;vdw-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
;bonded-lambda            = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

;fep-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
;coul-lambda              = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
;restraint-lambda         = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
;vdw-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
;bonded-lambda            = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

sc-alpha                 = 0.5 ; soft-core factor
"""

def compose_blocks(blocknames):
    """Compose desired blocks into a mdp file.

    ARGUMENTS
      blocknames (list of strings) - list of names of blocks to be concatenated

    RETURNS
      mdplines (list of strings) - concatenated blocks forming the contents of a gromacs .mdp file
    """

    contents = ""

    for blockname in blocknames:
        contents += mdpblock[blockname]

    mdplines = contents.split('\n')

    return mdplines
        
#=============================================================================================
# SUBROUTINES
#=============================================================================================
def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines
#=============================================================================================
def setupHydrationCalculation(solute, nreplicates = 1, verbose = True, jobname = "hydration"):
    """Set up an absolute alchemical hydration free energy calculation for the given molecule.

    ARGUMENTS
      solute (OEMol) - the molecule for which hydration free energy is to be computed (with fully explicit hydrogens) in the desired protonation state.

    OPTIONAL ARGUMENTS
      nreplicates (integer) - the number of replicates to set up (default: 1)
      verbose (boolean) - if True, extra debug information will be printed (default: True)
      jobname (string) - string to use for job name (default: "hydration")

    NOTES
      A directory will be created 'molecules/[molecule name]' as obtained from molecule.GetTitle().
    
    """

    # get current directory
    current_path = os.getcwd()

    # Center the solute molecule.
    OECenter(solute)

    # get molecule name
    solute_name = molecule.GetTitle()
    if verbose: print solute_name

    # create molecule path/directory
    work_path = os.path.abspath(os.path.join('molecules', solute_name))
    os.makedirs(work_path)

#    # Write mol2 file for the molecule.
#    writeMolecule(solute, os.path.join(solute_path, 'solute.mol2'))

#    # Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
#    parameterizeForGromacs(molecule, topology_filename = os.path.join(molecule_path,'solute.top'), coordinate_filename = os.path.join(molecule_path,'solute.gro'), charge_model = 'bcc', resname = 'MOL')

#    # Modify gromacs topology file for alchemical free energy calculation.
#    perturbGromacsTopology(os.path.join(molecule_path,'solute.top'), molecule)

#    # Convert gromacs topology to an .itp file suitable for inclusion.
#    top_to_itp(os.path.join(molecule_path,'solute.top'), os.path.join(molecule_path,'solute.itp'), moleculetype = molecule_name)
    
    # get pathnames
    mdrun = globals()['mdrun']
    grompp = globals()['grompp']
    editconf = globals()['editconf']    

    # SET UP SOLUTE TOPOLOGY

    if verbose: print "\nCONSTRUCTING SOLUTE TOPOLOGY"

    # Get formal charge of ligand.
    solute_charge = formalCharge(solute)
    if verbose: print "solute formal charge is %d" % solute_charge
    
    # Write molecule with explicit hydrogens to mol2 file.
    print "Writing solute mol2 file..."
    solute_mol2_filename = os.path.abspath(os.path.join(work_path, 'solute.mol2'))
    writeMolecule(solute, solute_mol2_filename)

    # Set substructure name (which will become residue name).
    print "Modifying molecule name..."
    modifySubstructureName(solute_mol2_filename, 'MOL')

    # Run antechamber to assign GAFF atom types.
    print "Running antechamber..."
    os.chdir(work_path)
    gaff_mol2_filename = os.path.join(work_path, 'solute.gaff.mol2')
    charge_model = 'bcc'
    command = 'antechamber -i %(solute_mol2_filename)s -fi mol2 -o solute.gaff.mol2 -fo mol2 -c %(charge_model)s -nc %(solute_charge)d > antechamber.out' % vars()
    if verbose: print command
    output = commands.getoutput(command)
    if verbose: print output
    os.chdir(current_path)

    # Generate frcmod file for additional GAFF parameters.
    solute_frcmod_filename = os.path.join(work_path, 'frcmod.solute')
    command = 'parmchk -i %(gaff_mol2_filename)s -f mol2 -o %(solute_frcmod_filename)s' % vars()
    if verbose: print command
    output = commands.getoutput(command)
    if verbose: print output

    # Run LEaP to generate topology / coordinates.
    solute_prmtop_filename = os.path.join(work_path,'solute.prmtop')
    solute_crd_filename = os.path.join(work_path,'solute.crd')
    solute_off_filename = os.path.join(work_path, 'solute.off')
    
    tleap_input_filename = os.path.join(work_path, 'setup-solute.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-solute.leap.out')
    contents = """
# Load GAFF parameters.
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# load solute
solute = loadMol2 %(gaff_mol2_filename)s

# check the solute
check solute

# report net charge
charge solute

# save AMBER parameters
saveAmberParm solute %(solute_prmtop_filename)s %(solute_crd_filename)s

# write .off file
saveOff solute %(solute_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    solute_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    solute_charge = int(round(float(solute_charge))) # round to nearest whole charge
    if verbose: print "solute charge is %d" % solute_charge    

    # PREPARE SOLVATED SOLUTE

    print "\nPREPARING SOLVATED SOLUTE"

    # create the directory if it doesn't exist
    solvent_path = os.path.join(work_path, 'solvent')
    if not os.path.exists(solvent_path):
        os.makedirs(solvent_path)

    # solvate the solute
    print "Solvating the solute with tleap..."
    system_prmtop_filename = os.path.join(solvent_path,'system.prmtop')
    system_crd_filename = os.path.join(solvent_path,'system.crd')
    tleap_input_filename = os.path.join(solvent_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(solvent_path, 'setup-system.leap.out')
    clearance = globals()['clearance'] # clearance around solute (in A)
    contents = """
source leaprc.ff99
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# Load solute.
loadOff %(solute_off_filename)s

# Create system.
system = combine { solute }
""" % vars()
    # add counterions
    if (solute_charge != 0):
        nions = abs(solute_charge)
        if solute_charge < 0: iontype = 'Na+'
        if solute_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
solvateBox system TIP3PBOX %(clearance)f iso

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()    
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)    

    # convert to gromacs
    print "Converting to gromacs..."
    amb2gmx = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'amb2gmx.pl')
    system_prefix = os.path.join(solvent_path, 'system')
    os.chdir(solvent_path)
    command = '%(amb2gmx)s --prmtop system.prmtop --crd system.crd --outname system' % vars()
    print command
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # Extract box size.
    g96_filename = os.path.join(solvent_path, 'system.g96')
    g96_file = open(g96_filename, 'r')
    g96_lines = g96_file.readlines()
    g96_file.close()
    box_size = zeros([3], float32)
    for line_number in range(len(g96_lines)):
        if g96_lines[line_number][0:3] == 'BOX':
            # parse line with box size in nm            
            line = g96_lines[line_number + 1]
            elements = line.split()
            box_size[0] = float(elements[0]) * 10.0
            box_size[1] = float(elements[1]) * 10.0
            box_size[2] = float(elements[2]) * 10.0            
    print "box_size = "
    print box_size

    # make a PDB file for checking
    print "Converting system to PDB..."
    os.chdir(solvent_path)
    command = 'cat system.crd | ambpdb -p system.prmtop > system.pdb' % vars()
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)
    
    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(solvent_path, 'system.top')
    lines = read_file(system_top_filename)
    perturb_atom_indices = list()
    indices = extract_section(lines, 'atoms')
    for index in indices:
        # extract the line
        line = stripcomments(lines[index])
        # parse the line
        elements = line.split()
        nelements = len(elements)
        # skip if not all elements found
        if (nelements < 8): continue
        # parse line
        atom = dict()
        atom['nr'] = int(elements[0])
        atom['type'] = elements[1]
        atom['resnr'] = int(elements[2])
        atom['residue'] = elements[3]
        atom['atom'] = elements[4]
        atom['cgnr'] = int(elements[5])
        atom['charge'] = float(elements[6])
        atom['mass'] = float(elements[7])
        # add those atoms in the solute to our list
        if atom['residue'] == 'MOL':
            perturb_atom_indices.append(atom['nr'])    
    perturbGromacsTopology(system_top_filename, solute, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = perturb_atom_indices)

    # set up replicates
    for replicate in range(nreplicates):
        # Create replicate directory.
        working_path = os.path.join(solvent_path, '%d' % replicate)
        os.makedirs(working_path)        

        # TODO: Modify solute coordinates in system.g96 to cycle through available conformations.
        
        # set up mdp files
        print "Writing mdp files for replicate %d..." % replicate

        mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'output']))
        mdpfile.write(os.path.join(working_path, 'minimize.mdp'))
        
        mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'output']))
        mdpfile.setParameter('nsteps', '10000') # 20 ps equilibration
        mdpfile.randomizeSeed() # randomize velocities
        mdpfile.write(os.path.join(working_path, 'equilibration.mdp'))

        mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'free-energy', 'output']))
        mdpfile.randomizeSeed() # randomize velocities
        mdpfile.write(os.path.join(working_path, 'production.mdp'))

        # write run script
        print "Writing run script..."
        solvent_path = os.path.abspath(solvent_path)
        contents = """\
set#!/bin/tcsh
#BSUB -J %(jobname)s_solvent
#BSUB -o %(working_path)s/outfile.%%J
#BSUB -e %(working_path)s/errfile.%%J
#BSUB -n 1
#BSUB -M 2000000
#BSUB -W 168:0

source $GMXRC
cd %(working_path)s

# constrained minimize
%(grompp)s -f minimize.mdp -c ../system.g96 -p ../system.top -o minimize.tpr -maxwarn 10000 -n ../system.ndx
%(mdrun)s -s minimize.tpr -x minimize.xtc -c minimize.g96 -e minimize.edr -g minimize.log

# equilibration
%(grompp)s -f equilibration.mdp -c minimize.g96 -p ../system.top -o equilibration.tpr -maxwarn 10000 -n ../system.ndx
%(mdrun)s -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.g96 -e equilibration.edr -g equilibration.log

# production
%(grompp)s -f production.mdp -c equilibration.g96 -p ../system.top -o production.tpr -maxwarn 10000 -n ../system.ndx
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.g96 -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
""" % vars()
        write_file(os.path.join(working_path, 'run.sh'), contents)

    # SET UP VACUUM SIMULATION

    # construct pathname for vacuum simulations
    vacuum_path = os.path.join(work_path, 'vacuum')
    if not os.path.exists(vacuum_path):
        os.makedirs(vacuum_path)

    # solvate the solute
    print "Preparing vacuum solute with tleap..."
    system_prmtop_filename = os.path.join(vacuum_path,'system.prmtop')
    system_crd_filename = os.path.join(vacuum_path,'system.crd')
    tleap_input_filename = os.path.join(vacuum_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(vacuum_path, 'setup-system.leap.out')
    clearance = 50.0 # clearance in A
    contents = """
source leaprc.ff99
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# Load solute.
loadOff %(solute_off_filename)s

# Create system.
system = combine { solute }
""" % vars()
    # add counterions
    if (solute_charge != 0):
        nions = abs(solute_charge)
        if solute_charge < 0: iontype = 'Na+'
        if solute_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """

# Create big box.
setBox system centers %(clearance)f iso

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()    
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)    
    
    # convert to gromacs
    print "Converting to gromacs..."    
    amb2gmx = os.path.join(os.getenv('MMTOOLSPATH'), 'converters', 'amb2gmx.pl')
    os.chdir(vacuum_path)    
    command = '%(amb2gmx)s --prmtop system.prmtop --crd system.crd --outname system' % vars()
    print command
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)    

    # make a PDB file for checking
    print "Converting system to PDB..."
    os.chdir(vacuum_path)
    command = 'cat system.crd | ambpdb -p system.prmtop > system.pdb' % vars()
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # write enlarged box for solvent because LEaP doesn't do it right.
    system_g96_filename = os.path.join(vacuum_path, 'system.g96')
    command = '%s -f %s -o %s -bt cubic -box %f %f %f -center 0.0 0.0 0.0' % (editconf, system_g96_filename, system_g96_filename, box_size[0]/10.0, box_size[1]/10.0, box_size[2]/10.0)
    print command
    output = commands.getoutput(command)
    print output

    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(vacuum_path, 'system.top')
    lines = read_file(system_top_filename)
    perturb_atom_indices = list()
    indices = extract_section(lines, 'atoms')
    nions = 0
    for index in indices:
        # extract the line
        line = stripcomments(lines[index])
        # parse the line
        elements = line.split()
        nelements = len(elements)
        # skip if not all elements found
        if (nelements < 8): continue
        # parse line
        atom = dict()
        atom['nr'] = int(elements[0])
        atom['type'] = elements[1]
        atom['resnr'] = int(elements[2])
        atom['residue'] = elements[3]
        atom['atom'] = elements[4]
        atom['cgnr'] = int(elements[5])
        atom['charge'] = float(elements[6])
        atom['mass'] = float(elements[7])
        # DEBUG
        # add those atoms in the solute to our list
        if atom['residue'] == 'MOL':
            perturb_atom_indices.append(atom['nr'])
        # add counterions to balance solute
        if ((atom['residue'] == 'Na+') or (atom['residue'] == 'Cl-')) and (nions < abs(solute_charge)):
            print "WARNING: Solute has net charge of %(solute_charge)d -- will perturb counterions too." % vars()
            perturb_atom_indices.append(atom['nr'])
            nions += 1
    perturbGromacsTopology(system_top_filename, solute, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = perturb_atom_indices)

    # construct replicates
    for replicate in range(nreplicates):
        # Create replicate path.
        working_path = os.path.join(vacuum_path, '%d' % replicate)
        os.makedirs(working_path)
    
        # TODO: Modify solute coordinates in system.g96 to cycle through available conformations.

        # construct .mdp files
        print "Writing mdp files for replicate %d..." % replicate

        mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-vacuum', 'constraints', 'output']))
        cutoff = box_size.min() / 10.0 / 2.0 - 0.0001 # compute cutoff
        mdpfile.setParameter('rlist', '%f' % cutoff)
        mdpfile.setParameter('rvdw', '%f' % cutoff)
        mdpfile.setParameter('rcoulomb', '%f' % cutoff)        
        mdpfile.write(os.path.join(working_path, 'minimize.mdp'))
        
        mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-vacuum', 'constraints', 'thermostat', 'free-energy', 'output']))
        mdpfile.randomizeSeed() # randomize velocities
        mdpfile.setParameter('rlist', '%f' % cutoff)
        mdpfile.setParameter('rvdw', '%f' % cutoff)
        mdpfile.setParameter('rcoulomb', '%f' % cutoff)        
        mdpfile.write(os.path.join(working_path, 'production.mdp'))

        # write shell script for minimization and production
        contents = """\
#!/bin/tcsh
#BSUB -J %(jobname)s_vacuum
#BSUB -o %(working_path)s/outfile.%%J
#BSUB -e %(working_path)s/errfile.%%J
#BSUB -n 1
#BSUB -M 2000000
#BSUB -W 168:0

source $GMXRC
cd %(working_path)s

# constrained minimize
%(grompp)s -f minimize.mdp -c ../system.g96 -p ../system.top -o minimize.tpr -maxwarn 10000 -n ../system.ndx
%(mdrun)s -s minimize.tpr -x minimize.xtc -c minimize.g96 -e minimize.edr -g minimize.log

# production
%(grompp)s -f production.mdp -c minimize.g96 -p ../system.top -o production.tpr -maxwarn 10000 -n ../system.ndx
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.g96 -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
""" % vars()    
        write_file(os.path.join(working_path, 'run.sh'), contents)

    return


#=============================================================================================
# MAIN
#=============================================================================================

# Create a molecule.
# molecule = createMoleculeFromIUPAC('phenol')
# molecule = createMoleculeFromIUPAC('biphenyl')
# molecule = createMoleculeFromIUPAC('4-(2-Methoxyphenoxy) benzenesulfonyl chloride')

# molecule = readMolecule('/Users/jchodera/projects/CUP/SAMPL-2008/jnk3/data-from-openeye/jnk.aff/jnk.aff-1.sdf', normalize = True)

# List of tuples describing molecules to construct:
# ('common name to use for directory', 'IUPAC name', formal_charge to select)
molecules = [
#    ('phenol', 'phenol', 0), # phenol
#    ('N-methylacetamide', 'N-methylacetamide', 0),
#    ('E-oct-2-enal', 'E-oct-2-enal', 0),
#    ('catechol', 'catechol', 0),
#    ('trimethylphosphate', 'trimethylphosphate', 0),
#    ('trimethylphosphate', 'triacetylglycerol', 0),
#    ('imatinib', '4-[(4-methylpiperazin-1-yl)methyl]-N-[4-methyl-3-[(4-pyridin-3-ylpyrimidin-2-yl)amino]-phenyl]-benzamide', 0), # imatinib
#    ('NanoKid', '2-(2,5-bis(3,3-dimethylbut-1-ynyl)-4-(2-(3,5-di(pent-1-ynyl)phenyl)ethynyl)phenyl)-1,3-dioxolane', 0), # NanoKid - http://en.wikipedia.org/wiki/Nanoputian
#    ('lipitor', '[R-(R*,R*)]-2-(4-fluorophenyl)-beta,delta-dihydroxy-5-(1-methylethyl)-3-phenyl-4-[(phenylamino)carbonyl]-1H-pyrrole-1-heptanoic acid', 0), # broken

    # trypsin inhibitors from Alan E. Mark study JACS 125:10570, 2003.
    ('benzamidine', 'benzamidine', 1), # benzamidine
    ('p-methylbenzamidine', 'p-methylbenzamidine', 1), # R = Me
    ('p-ethylbenzamidine', 'p-ethylbenzamidine', 1), # R = Et
    ('p-n-propylbenzamidine', 'p-n-propylbenzamidine', 1), # R = n-Pr
    ('p-isopropylbenzamidine', 'p-isopropylbenzamidine', 1), # R = i-Pr
    ('p-n-butylbenzamidine', 'p-n-butylbenzamidine', 1), # R = n-Bu
    ('p-t-butylbenzamidine', 'p-t-butylbenzamidine', 1), # R = t-Bu
    ('p-n-pentylbenzamidine', 'p-n-pentylbenzamidine', 1), # R = n-Pent
    ('p-n-hexylbenzamidine', 'p-n-hexylbenzamidine', 1), # R = n-Hex    
    ]

for (molecule_name, molecule_iupac_name, formal_charge) in molecules:
    print "Setting up %s" % molecule_name

    # Create the molecule from the common name with the desired formal charge.
    molecule = createMoleculeFromIUPAC(molecule_iupac_name, charge = formal_charge)
    if not molecule:
        print "Could not build '%(molecule_name)s' from IUPAC name '%(molecule_iupac_name)s', skipping..." % vars()
        continue
    print "Net charge is %d" % formalCharge(molecule)

    # Replace the title with the common name
    molecule.SetTitle(molecule_name)

    # Expand set of conformations so we have multiple conformations to start from.
    expandConformations(molecule)
    
    # Set up multiple replicates of the hydration free energy calculation.
    setupHydrationCalculation(molecule, nreplicates = nreplicates, jobname = molecule_name)

    print "\n\n"

