#=============================================================================================
# Example illustrating use of 'moltools' in computing hydration free energies of small molecules.
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

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

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
nsteps                   = 1000000 ; 2 ns
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
nstlog                   = 20 ; (must be same as nstdgdl for analysis scripts)
nstenergy                = 20 ; (must be same as nstdgdl for analysis scripts)
; Output frequency and precision for xtc file = 
nstxtcout                = 20 ; change this to larger value later
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
vdw-type                 = switch
; cut-off lengths        = 
rvdw-switch              = 0.8
rvdw                     = 0.9
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
nstfep                   = 20 ; 0.1 ps between weight updates (must be integer multiple of nstlist)
nstdgdl                  = 20 ; 0.1 ps between writing energies (must be same as nstdgdl for analysis scripts)

; weight update scheme
; lambda-mc                = gibbs-wang-landau ; Wang-Landau with waste recycling
lambda-mc                = barker-transition ; UBAR
mc-wldelta               = 1.0 ; initial delta factor for Wang-Landau (in kT)
mc-wlscale               = 0.5 ; scalar by which delta is scaled for Wang-Landau
mc-nratio                = 0.8 ; flatness criterion -- histograms are reset after all states are sampled within mc-nratio factor of the mean

; state transition probability
move-mc                  = metropolized-gibbs ; Metropolized Gibbs for fastest mixing of states

; starting and stopping
mc-nstart                = 100 ; number of updates to perform per state for driving through each state
mc-nequil                = 100000 ; number of steps before freezing weights (200 ps)

init-lambda              = 1 ; initial state

; schedule for switching off lambdas
; first, restraints are turned on as charges are switched off
; next, vdw and torsions are switched off
fep-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.0 0.00 0.0 0.00 0.0 0.0 ; for global scaling (don't need)
coul-lambda              = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling electrostatics
restraint-lambda         = 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.00 1.0 1.00 1.0 1.00 1.0 1.0 ; for scaling restraints
vdw-lambda               = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling vdw interactions
bonded-lambda            = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 ; for scaling torsions

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
#=============================================================================================
def setupHydrationCalculation(molecule):
    """
    """

    # get current directory
    current_path = os.getcwd()

    # Center the molecule.
    OECenter(molecule)

    # get molecule name
    molecule_name = molecule.GetTitle()
    print molecule_name

    # create molecule path/directory
    molecule_path = os.path.abspath(os.path.join('molecules', molecule_name))
    os.makedirs(molecule_path)

    # Write mol2 file for the molecule.
    writeMolecule(molecule, os.path.join(molecule_path, 'solute.mol2'))

    # Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
    parameterizeForGromacs(molecule, topology_filename = os.path.join(molecule_path,'solute.top'), coordinate_filename = os.path.join(molecule_path,'solute.gro'), charge_model = 'bcc', resname = 'MOL')

    # Modify gromacs topology file for alchemical free energy calculation.
    perturbGromacsTopology(os.path.join(molecule_path,'solute.top'), molecule)

    # Convert gromacs topology to an .itp file suitable for inclusion.
    top_to_itp(os.path.join(molecule_path,'solute.top'), os.path.join(molecule_path,'solute.itp'), moleculetype = molecule_name)

    # SET UP VACUUM SIMULATION

    # construct pathname for vacuum simulations
    vacuum_path = os.path.join(molecule_path, 'vacuum')
    os.mkdir(vacuum_path)
    os.chdir(vacuum_path)
    
    # create molecule in enlarged box
    command = 'editconf_d -f ../solute.gro -o system.gro -bt cubic -d 50.0 -center 0.0 0.0 0.0'
    print command
    output = commands.getoutput(command)
    print output

    # create topology file
    topology_contents = """
; amber forcefield
#include "ffamber03.itp"

; my molecule
#include "../solute.itp"

[ system ]
%(molecule_name)s in vacuum

[ molecules ]
; Compound        nmols
%(molecule_name)s            1
""" % vars()
    writeFile('system.top', topology_contents)

    # construct .mdp files
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-vacuum', 'output']))
    mdpfile.write(os.path.join(vacuum_path, 'minimization-unconstrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-vacuum', 'constraints', 'output']))
    mdpfile.write(os.path.join(vacuum_path, 'minimization-constrained.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-vacuum', 'constraints', 'thermostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(vacuum_path, 'production.mdp'))

    # write shell script for minimization and production
    contents = """\
#/bin/tcsh

source /Users/jchodera/local/gromacs-3.1.4-fe/i386-apple-darwin8.11.1/bin/GMXRC.csh
setenv GROMPP grompp_d
setenv MDRUN mdrun_d
    
# unconstrained minimization
$GROMPP -f minimization-unconstrained.mdp -c system.gro -p system.top -o minimization-unconstrained.tpr -maxwarn 10000
$MDRUN -s minimization-unconstrained.tpr -x minimization-unconstrained.xtc -c minimized-unconstrained.gro -e minimization-unconstrained.edr -g minimization-unconstrained.log
    
# constrained minimization
$GROMPP -f minimization-constrained.mdp -c minimized-unconstrained.gro -p system.top -o minimization-constrained.tpr -maxwarn 10000
$MDRUN -s minimization-constrained.tpr -x minimization-constrained.xtc -c minimized-constrained.gro -e minimization-constrained.edr -g minimization-constrained.log

# production
$GROMPP -f production.mdp -c minimized-constrained.gro -p system.top -o production.tpr -maxwarn 10000
$MDRUN -v -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log
echo 0 | trjconv_d -f production.xtc -o production.pdb -s production.tpr
""" % vars()    
    writeFile(os.path.join(vacuum_path, 'run.sh'), contents)

    # SET UP SOLVATED CALCULATION

    # construct pathname for solvent simulations
    solvated_path = os.path.join(molecule_path, 'solvent')
    os.mkdir(solvated_path)
    os.chdir(solvated_path)
    
    # write enlarged box for solvent
    command = 'editconf_d -f ../solute.gro -o solute.gro -bt octahedron -d 1.0 -center 0.0 0.0 0.0'
    print command
    output = commands.getoutput(command)
    print output

    # write topology
    topology_contents = """
; amber forcefield
#include "ffamber03.itp"

; my molecule
#include "../solute.itp"

; TIP3P water
#include "ffamber_tip3p.itp"

[ system ]
%(molecule_name)s in solvent

[ molecules ]
; Compound        nmols
%(molecule_name)s            1
""" % vars()
    writeFile('system.top', topology_contents)

    # solvate
    command = 'genbox_d -cp solute.gro -cs ffamber_tip3p.gro -o system.gro -p system.top'
    # command = 'genbox_d -cp solute.gro -cs tip3p-pme-waterbox.gro -o system.gro -p system.top' # use our special tip3p box
    # command = 'genbox_d -cp solute.gro -cs tip3p-pme-waterbox.gro -o system.gro -p system.top -vdwd 0.0' # insert waters but don't cull overlap
    print command
    output = commands.getoutput(command)
    print output
    
    # construct .mdp files
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'output']))
    mdpfile.write(os.path.join(solvated_path, 'minimization-unconstrained.mdp'))
    
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'output']))
    mdpfile.write(os.path.join(solvated_path, 'minimization-constrained.mdp'))
    
#    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat-equilibration', 'barostat', 'output']))
#    mdpfile.randomizeSeed() # randomize velocities
#    mdpfile.setParameter('nsteps', '25000') # 10 ps equilibration
#    mdpfile.write(os.path.join(solvated_path, 'equilibration.mdp'))
    
    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(solvated_path, 'production.mdp'))

    # write run script
    contents = """\
#/bin/tcsh

source /Users/jchodera/local/gromacs-3.1.4-fe/i386-apple-darwin8.11.1/bin/GMXRC.csh
setenv GROMPP grompp_d
setenv MDRUN mdrun_d
    
# unconstrained minimization
$GROMPP -f minimization-unconstrained.mdp -c system.gro -p system.top -o minimization-unconstrained.tpr -maxwarn 10000
$MDRUN -s minimization-unconstrained.tpr -x minimization-unconstrained.xtc -c minimized-unconstrained.gro -e minimization-unconstrained.edr -g minimization-unconstrained.log

# constrained minimization
$GROMPP -f minimization-constrained.mdp -c minimized-unconstrained.gro -p system.top -o minimization-constrained.tpr -maxwarn 10000
$MDRUN -s minimization-constrained.tpr -x minimization-constrained.xtc -c minimized-constrained.gro -e minimization-constrained.edr -g minimization-constrained.log

# equilibration
#$GROMPP -f equilibration.mdp -c minimized-constrained.gro -p system.top -o equilibration.tpr -maxwarn 10000
#$MDRUN -v -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.gro -e equilibration.edr -g equilibration.log
#echo 0 | trjconv_d -f equilibration.xtc -o equilibration.pdb -s equilibration.tpr

# production
$GROMPP -f production.mdp -c minimized-constrained.gro -p system.top -o production.tpr -maxwarn 10000
#$GROMPP -f production.mdp -c equilibration.gro -p system.top -o production.tpr -maxwarn 10000
#$GROMPP -f production.mdp -c system.gro -p system.top -o production.tpr -maxwarn 10000
$MDRUN -v -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log
echo 0 | trjconv_d -f production.xtc -o production.pdb -s production.tpr
""" % vars()
    writeFile(os.path.join(solvated_path, 'run.sh'), contents)

    # restore working directory
    os.chdir(current_path)

    return


#=============================================================================================
# MAIN
#=============================================================================================

# Create a molecule.
# molecule = createMoleculeFromIUPAC('phenol')
# molecule = createMoleculeFromIUPAC('biphenyl')
# molecule = createMoleculeFromIUPAC('4-(2-Methoxyphenoxy) benzenesulfonyl chloride')

# molecule = readMolecule('/Users/jchodera/projects/CUP/SAMPL-2008/jnk3/data-from-openeye/jnk.aff/jnk.aff-1.sdf', normalize = True)

molecules = [
#    ('phenol', 'phenol'), # phenol
#    ('N-methylacetamide', 'N-methylacetamide'),
#    ('E-oct-2-enal', 'E-oct-2-enal'),
#    ('catechol', 'catechol'),
#    ('trimethylphosphate', 'trimethylphosphate'),
#    ('trimethylphosphate', 'triacetylglycerol'),
#    ('imatinib', '4-[(4-methylpiperazin-1-yl)methyl]-N-[4-methyl-3-[(4-pyridin-3-ylpyrimidin-2-yl)amino]-phenyl]-benzamide'), # imatinib
#    ('NanoKid', '2-(2,5-bis(3,3-dimethylbut-1-ynyl)-4-(2-(3,5-di(pent-1-ynyl)phenyl)ethynyl)phenyl)-1,3-dioxolane'), # NanoKid - http://en.wikipedia.org/wiki/Nanoputian
#    ('lipitor', '[R-(R*,R*)]-2-(4-fluorophenyl)-beta,delta-dihydroxy-5-(1-methylethyl)-3-phenyl-4-[(phenylamino)carbonyl]-1H-pyrrole-1-heptanoic acid'), # broken
    ('benzamidine', 'benzamidine'), # benzamidine
    ]

for (molecule_name, molecule_iupac_name) in molecules:
    print molecule_name

    # Create the molecule from the common name
    molecule = createMoleculeFromIUPAC(molecule_iupac_name)

    # Replace the title with the common name
    molecule.SetTitle(molecule_name)

    # Set up hydration free energy calculation
    setupHydrationCalculation(molecule)

    print "\n\n"

