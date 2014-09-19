#=============================================================================================
# Example illustrating use of 'moltools' to set up free energy calculation of annihilating a ligand in a protein-ligand complex.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-01-23
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.moltools.ligandtools import *
from mmtools.gromacstools import *
from mmtools.gromacstools.MdpFile import *
from mmtools.gromacstools.System import *
import commands
import os    

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

nreplicates_solvent = 1 # number of replicates to set up of the solvated ligand
nreplicates_complex = 1 # number of replicates to set up of the complexed ligand

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
nsteps                   = 100000 ; 200 ps
; mode for center of mass motion removal = 
comm-mode                = Linear
; number of steps for center of mass motion removal = 
nstcomm                  = 500
; group(s) for center of mass motion removal = 
comm-grps                = 
""" 

# output control block
mdpblock['output'] = """
; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 500 ; (must be integral multiple of nstdgdl for post-analysis)
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 50 ; (must be same as nstdgdl for analysis scripts)
nstenergy                = 50 ; (must be same as nstdgdl for analysis scripts)
; Output frequency and precision for xtc file = 
nstxtcout                = 500 ; change this to larger value later
xtc-precision            = 1000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
xtc-grps                 = System
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
gen_temp                 = 298.0 ; K
gen_seed                 = 12345

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Temperature coupling   = 
Tcoupl                   = AndersenII
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tau_t                    = 5.0
ref_t                    = 298.0 ; K
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

mdpblock['nonbonded-solvent'] = """
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

# free energy block
mdpblock['free-energy'] = """
; OPTIONS FOR EXPANDED ENSEMBLE SIMULATIONS
; Free energy control stuff = 
free-energy              = mutate ; annihilate electrostatics and Lennard-Jones
nstfep                   = 50 ; 0.1 ps between weight updates (must be integer multiple of nstlist)
nstdgdl                  = 50 ; 0.1 ps between writing energies (must be same as nstdgdl for analysis scripts)

; weight update scheme
lambda-mc                = gibbs-wang-landau ; Wang-Landau with waste recycling
mc-wldelta               = 1.0 ; initial delta factor for Wang-Landau (in kJ/mol)
mc-wlscale               = 0.5 ; scalar by which delta is scaled for Wang-Landau
mc-nratio                = 0.8 ; flatness criterion (between 0 and 1)

; state transition probability
move-mc                  = metropolized-gibbs ; Metropolized Gibbs for fastest mixing of states

; starting and stopping
mc-nstart                = 0 ; number of updates to perform per state for driving through each state
mc-nequil                = 1000 ; number of updates before freezing weights (50 ps)

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

def readFile(filename):
    """ Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines


def setup_binding_free_energy(protein_pdb_filename, ligand_filename, work_path, nreplicates_solvent = 1, nreplicates_complex = 1):
    """Set up generalized ensemble simulation to calculate a protein-ligand binding free energy.

    ARGUMENTS
      protein_pdb_filename (string) - filename of PDB file for the protein, without any missing heavy atoms and with appropriate protonation state already assigned
      ligand_filename (string) - filename of some type OpenEye can read (e.g. mol2) describing the ligand with suitable docked coordinates
      work_path (string) - directory to place files for binding free energy calculation (created if does not already exist)

    OPTIONAL ARGUMENTS
      nreplicates_solvent (int) - number of replicates to set up of solvated ligand (default: 1)
      nreplicates_complex (int) - number of replicates to set up of the complexed ligand (default: 1)
    """

    # parameters
    forcefield = 'ffamber96' # forcefield
    box_clearance_in_nm = 0.3 # clearance from system to box edges (in nm) # DEBUG set to 1.0

    # gromacs filenames
    pdb2gmx = 'pdb2gmx_d'
    editconf = 'editconf_d'
    genbox = 'genbox_d'
    grompp = 'grompp_d'
    mdrun = 'mdrun_d'

    # Create directory for binding free energy calculation if it doesn't exist.
    #DEBUG os.makedirs(work_path)

    # SET UP PROTEIN TOPOLOGY

    print "\nConstructing protein topology..."

    # Create .top and .gro files for protein using pdb2gmx and desired forcefield
    # look up forcefield index
    lines = readFile(os.path.join(os.environ['GMXDATA'], 'top', 'FF.dat'))
    for forcefield_index in range(1, len(lines)):
        if lines[forcefield_index].find(forcefield) != -1: break
    forcefield_index -= 1
    # run pdb2gmx
    protein_gro_filename = os.path.join(work_path, 'protein.gro')
    protein_top_filename = os.path.join(work_path, 'protein.top')
    command = 'echo %(forcefield_index)d | %(pdb2gmx)s -f %(protein_pdb_filename)s -o %(protein_gro_filename)s -p %(protein_top_filename)s -ignh -H14' % vars()
    print command
    #DEBUG output = commands.getoutput(command)
    #DEBUG print output

    # SET UP LIGAND TOPOLOGY

    print "\nConstructing ligand topology..."

    # Read the ligand into OEMol.
    ligand = readMolecule(ligand_filename, normalize = True)
    
    # Generate GAFF parameters for gromacs for small molecule, using antechamber to generate AM1-BCC charges.
    ligand_top_filename = os.path.join(work_path,'ligand.top')
    ligand_gro_filename = os.path.join(work_path,'ligand.gro')
    #DEBUG parameterizeForGromacs(ligand, topology_filename = ligand_top_filename, coordinate_filename = ligand_gro_filename, charge_model = 'bcc', resname = 'LIG')

    # Modify gromacs topology file for alchemical free energy calculation.
    #DEBUG perturbGromacsTopology(os.path.join(work_path,'ligand.top'), ligand)
    
    # SET UP SOLVATED LIGAND

    print "\nSetting up solvated ligand..."

    # Create directory to contain ligand in solution.
    solvated_ligand_path = os.path.join(work_path, 'solvated-ligand')
    #DEBUG os.makedirs(solvated_ligand_path)

    # Convert gromacs topology to an .itp file suitable for inclusion.
    top_to_itp(os.path.join(work_path,'ligand.top'), os.path.join(solvated_ligand_path,'solute.itp'), moleculetype = 'solute')

    # write enlarged box for solvent
    system_gro_filename = os.path.join(solvated_ligand_path, 'system.gro')
    system_top_filename = os.path.join(solvated_ligand_path, 'system.top')
    command = '%(editconf)s -f %(ligand_gro_filename)s -o %(system_gro_filename)s -bt octahedron -d %(box_clearance_in_nm)f' % vars()
    print command
    output = commands.getoutput(command)
    print output

    # write topology for system
    contents = """
; amber forcefield
#include "ffamber03.itp"

; my molecule
#include "solute.itp"

; TIP3P water
#include "ffamber_tip3p.itp"

[ system ]
solute

[ molecules ]
; Compound        nmols
solute            1
""" % vars()
    writeFile(system_top_filename, contents)

    # solvate
    command = '%(genbox)s -cp %(system_gro_filename)s -cs ffamber_tip3p.gro -o %(system_gro_filename)s -p %(system_top_filename)s' % vars()
    print command
    output = commands.getoutput(command)
    print output
    
    # construct .mdp files
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'output']))
    mdpfile.write(os.path.join(solvated_ligand_path, 'minimization-unconstrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'output']))
    mdpfile.write(os.path.join(solvated_ligand_path, 'minimization-constrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(solvated_ligand_path, 'production.mdp'))

    # write run script
    contents = """\
#/bin/tcsh

source /Users/jchodera/local/gromacs-3.1.4-fe/i386-apple-darwin8.11.1/bin/GMXRC.csh
    
# unconstrained minimization
%(grompp)s -f minimization-unconstrained.mdp -c system.gro -p system.top -o minimization-unconstrained.tpr -maxwarn 10000
%(mdrun)s -s minimization-unconstrained.tpr -x minimization-unconstrained.xtc -c minimized-unconstrained.gro -e minimization-unconstrained.edr -g minimization-unconstrained.log

# constrained minimization
%(grompp)s -f minimization-constrained.mdp -c minimized-unconstrained.gro -p system.top -o minimization-constrained.tpr -maxwarn 10000
%(mdrun)s -s minimization-constrained.tpr -x minimization-constrained.xtc -c minimized-constrained.gro -e minimization-constrained.edr -g minimization-constrained.log

# production
%(grompp)s -f production.mdp -c minimized-unconstrained.gro -p system.top -o production.tpr -maxwarn 10000
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log
echo 0 | trjconv_d -f production.xtc -o production.pdb -s production.tpr
""" % vars()
    writeFile(os.path.join(solvated_ligand_path, 'run.sh'), contents)

    # SET UP SOLVATED COMPLEX

    # Create directory to contain complex.
    complex_path = os.path.join(work_path, 'solvated-complex')
    #DEBUG os.makedirs(complex_path)

    # Merge protein and liand .top files.
    complex_top_filename = os.path.join(complex_path, 'complex.top')
    complex_gro_filename = os.path.join(complex_path, 'complex.gro')
    merge_protein_ligand_topologies(protein_top_filename, protein_gro_filename, ligand_top_filename, ligand_gro_filename, complex_top_filename, complex_gro_filename)

    # Alter topology file to use ffamber TIP3P.
    # system.useTIP3P(complex_top_filename)

    # Determine total charge.
    total_charge = totalCharge(complex_top_filename)
    print "total_charge = %f" % total_charge
    
    # Create a System object to help us out.
    system = System(protein_pdb_filename, finalOutputDir = complex_path, finalOutputName = "protein", useff = forcefield)
    system.setup.setSaltConditions('MgCl2', 0.010)
    system.setup.set_boxType = 'octahedron'
    system.setup.set_boxSoluteDistance(box_clearance_in_nm)
    system.files.topfile = complex_top_filename
    system.totalChargeBeforeIons = total_charge

    # Create a periodic boundary conditions box.
    command = '%(editconf)s -bt octahedron -f %(complex_gro_filename)s -o %(complex_gro_filename)s -d %(box_clearance_in_nm)f' % vars()
    output = commands.getoutput(command)    
    print output

    # solvate the box with TIP3P water
    command = '%(genbox)s -cp %(complex_gro_filename)s -cs ffamber_tip3p.gro -o %(complex_gro_filename)s -p %(complex_top_filename)s' % vars()
    output = commands.getoutput(command)
    print output    

    # calculate how many ions to add to the box
    [np, nn, nwaters] = system.counterions()	 
    print "Adding %(np)d positive and %(nn)d negative to %(nwaters)d waters" % vars()
    

    # minimize the system
    ### make a tpr file with grompp
    #grompp = '%s/grompp -f %s -c %s -o %s -p %s '%(os.environ['GMXPATH'], self.files.mdpfile_Minimization, self.files.grofile, self.files.tprfile, self.files.topfile)
    #self.rungmx( grompp, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )	

    ### run minimization
    #minimize = '%s/mdrun -v -s %s -c %s '%(os.environ['GMXPATH'], self.files.tprfile, self.files.next_gro() )
    #self.rungmx( minimize, mockrun=self.mockrun, checkForFatalErrors=self.checkForFatalErrors )
    #self.files.increment_gro()    # must increment filename for any new gmx file 
    
    # do the rest of the preparation steps
    #self.postSolvationPreparationSteps()


    # construct .mdp files
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'output']))
    mdpfile.write(os.path.join(complex_path, 'minimization-unconstrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'output']))
    mdpfile.write(os.path.join(complex_path, 'minimization-constrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(complex_path, 'production.mdp'))    

    # write run script for minimization and production
    contents = """\
#/bin/tcsh

source /Users/jchodera/local/gromacs-3.1.4-fe/i386-apple-darwin8.11.1/bin/GMXRC.csh
    
# unconstrained minimization
%(grompp)s -f minimization-unconstrained.mdp -c system.gro -p system.top -o minimization-unconstrained.tpr -maxwarn 10000
%(mdrun)s -s minimization-unconstrained.tpr -x minimization-unconstrained.xtc -c minimized-unconstrained.gro -e minimization-unconstrained.edr -g minimization-unconstrained.log

# constrained minimization
%(grompp)s -f minimization-constrained.mdp -c minimized-unconstrained.gro -p system.top -o minimization-constrained.tpr -maxwarn 10000
%(mdrun)s -s minimization-constrained.tpr -x minimization-constrained.xtc -c minimized-constrained.gro -e minimization-constrained.edr -g minimization-constrained.log

# production
%(grompp)s -f production.mdp -c minimized-unconstrained.gro -p system.top -o production.tpr -maxwarn 10000
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log
echo 0 | trjconv_d -f production.xtc -o production.pdb -s production.tpr
""" % vars()
    writeFile(os.path.join(complex_path, 'run.sh'), contents)



    
         
    # Convert gromacs topology to an .itp file suitable for inclusion.
    # top_to_itp(os.path.join(molecule_path,'solute.top'), os.path.join(molecule_path,'solute.itp'), moleculetype = molecule_name)

    return


#=============================================================================================
# MAIN
#=============================================================================================

# set pathnames for testing
protein_pdb_filename = 'complex-structures/j3_1_automodel_mcce_out.pdb'
ligand_filename = 'complex-structures/jnk.aff-1.sdf'
work_path = 'jnk.aff-1'

# set up free energy calculation
setup_binding_free_energy(protein_pdb_filename, ligand_filename, work_path, nreplicates_solvent=nreplicates_solvent, nreplicates_complex=nreplicates_complex)

