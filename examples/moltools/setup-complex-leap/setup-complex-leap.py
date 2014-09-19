#=============================================================================================
# Set up a system for alchemical free energy calculations in gromacs using AMBER leap.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-01
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
import os.path
import shutil

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

# JDC LAPTOP
#parameter_path = "parameters/"
#protein_path = 'source-structures/unphosphorylated-cocrystal-structures' # source of cocrystal structures
#ligand_path = 'source-structures/inhibitors' # source of ligand structures
#GMXRC = '/Users/jchodera/local/gromacs-3.1.4-fe/i386-apple-darwin8.11.1/bin/GMXRC.csh'
#grompp = 'grompp_d'
#mdrun = 'mdrun_d'
#editconf = 'editconf_d'

# IMRAN CLUSTER
parameter_path = "/home/ihaque/buildtools/mmtools/examples/moltools/setup-complex-leap/parameters/" # path where additional parameters can be found
protein_path = '/home/ihaque/sampl/automodel_phos_fixed'
ligand_path = '/home/ihaque/sampl/srcfiles'
GMXRC = '/home/ihaque/local/gromacs/3.1.4-fe/bin/GMXRC'
grompp = 'grompp'
mdrun = 'mdrun'
editconf = 'editconf'

compute_angles = '/home/ihaque/buildtools/mmtools/examples/moltools/setup-complex-leap/computeangles.py'
restraint_topology = '/home/ihaque/buildtools/mmtools/examples/moltools/setup-complex-leap/restraint_topology.py'

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
nstlog                   = 20 ; 
nstenergy                = 20 ;
; Output frequency and precision for xtc file = 
nstxtcout                = 50 ;
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
    
mdpblock['restraints'] = """
;Enable dihedral and distance restraints
dihre=simple
dihre-fc=1
nstdihreout=1000
disre=simple
disre_fc=1
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
rvdw-switch              = 0.85 ; MRS has 0.7
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
optimize_fft             = yes
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
vdw-type                 = cut-off
; cut-off lengths        = 
rvdw                     = 23.0 ; to emulate no cutoff
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
mc-nstart                = 200 ; number of updates to perform per state for driving through each state
mc-nequil                = 250000 ; number of steps before freezing weights (200 ps)

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
def rename_pdb_for_amber(protein_pdb_filename, protein_amberpdb_filename):
    """Rename residues in a PDB file to match AMBER naming convention, omitting atoms AMBER will not recognize.

    ARGUMENTS
      protein_pdb_filename (string) - name of source PDB file to be read
      protein_amberpdb_filename (string) - name of PDB file to be written with AMBER naming conventions

    NOTE
      This conversion is still jnk3-specific at this time.
      Four-letter residue names aren't properly parsed, but a hack is in.

    """
    
    # read the PDB file
    lines = read_file(protein_pdb_filename)

    # count phosphate hydrogens
    count_tpo_hydrogens = int(commands.getoutput('grep -c HO.P.TPO %(protein_pdb_filename)s' % vars()))
    count_ptr_hydrogens = int(commands.getoutput('grep -c HO.P.PTR %(protein_pdb_filename)s' % vars()))

    # allocate storage for target file contents
    outlines = list()

    # parse PDB file, creating new contents
    for line in lines:
        # only ATOM records will be retained
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            serial = line[6:11]
            name = line[12:16]
            altLoc = line[16:17]
            resname = line[17:21]
            chainid = line[21:22]
            resseq = line[22:26]
            idcode = line[26:27]
            remainder = line[27:]

            # drop the line it is a hydrogen
            if name[1] == 'H': continue
            # phosphorylated residues have a strange naming convention on their hydrogens
            if name[0] == 'H': continue

            # make residue name substitutions
            if resname == 'LYP ': resname = 'LYS '
            if resname == 'CYN ': resname = 'CYS '
            if resname == 'NSER': resname = 'SER '
            if (resname == 'CGLU') or (resname == 'CGLH'):
                resname = 'GLU '
                # rename terminal oxygens
                if name == ' OC1': name = ' OXT'
                if name == ' OC2': name = ' O  '
            
            # determine phosphate naming
            if resname == 'TPO ':
                if name == ' OG1': name = ' OG '
                if count_tpo_hydrogens == 2: resname = 'T1P '
                else: resname = 'T2P '                    
            if resname == 'PTR ':
                if name == ' OH ': name = ' OG '
                if count_ptr_hydrogens == 2: resname = 'Y1P '
                else: resname = 'Y2P '

            # renumber serial
            serial = '% 5d' % (len(outlines) + 1)
            
            # reconstitute line
            outline = 'ATOM  ' + serial + ' ' + name + altLoc + resname + chainid + resseq + idcode + remainder
            outlines.append(outline)

    write_file(protein_amberpdb_filename, outlines)
            
    return
#=============================================================================================
def setup_system(protein_pdb_filename, ligand_filename, work_path, parameter_path, GMXRC, jobname):
    """Set up a system for alchemical free energy calculation in gromacs.

    ARGUMENTS
      protein_pdb_filename (string) - name of ffamber-named protein PDB file
      ligand_filename (string) - name of mol2 file describing ligand
      work_path (string) - name of directory to place files in

    """

    # get current directory
    current_path = os.getcwd()

    # create the work directory if it doesn't exist
    if not os.path.exists(work_path):
        os.makedirs(work_path)

    mdrun = globals()['mdrun']
    grompp = globals()['grompp']

    # SET UP PROTEIN TOPOLOGY

    print "\nCONSTRUCTING PROTEIN TOPOLOGY"

    # convert protein PDB file to AMBER naming conventions, dropping atoms that AMBER won't recognize (like protons)
    print "Converting PDB naming to AMBER..."
    protein_amberpdb_filename = os.path.join(work_path, 'protein-ambernames.pdb')
    rename_pdb_for_amber(protein_pdb_filename, protein_amberpdb_filename)

    # run leap to set up protein and report on net charge
    print "Running LEaP to set up protein..."

    protein_prmtop_filename = os.path.join(work_path,'protein.prmtop')
    protein_crd_filename = os.path.join(work_path,'protein.crd')
    protein_off_filename = os.path.join(work_path, 'protein.off')

    tleap_input_filename = os.path.join(work_path, 'setup-protein.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-protein.leap.out')

    contents = """
# Load AMBER ff99 parameters.
source leaprc.ff99

# Load phosphoresidue parameters.
loadoff %(parameter_path)s/T1P.off
loadoff %(parameter_path)s/T2P.off
loadoff %(parameter_path)s/Y1P.off
loadoff %(parameter_path)s/Y2P.off
loadamberparams %(parameter_path)s/frcmod_t1p
loadamberparams %(parameter_path)s/frcmod_t2p
loadamberparams %(parameter_path)s/frcmod_y1p
loadamberparams %(parameter_path)s/frcmod_y2p

# Load PDB file with AMBER naming conventions, stripped of hydrogens.
protein = loadpdb %(protein_amberpdb_filename)s

# check the system
check protein

# report net charge
charge protein

# write out parameters
saveAmberParm protein %(protein_prmtop_filename)s %(protein_crd_filename)s

# write as LEaP object
saveOff protein %(protein_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    protein_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    protein_charge = int(round(float(protein_charge))) # round to nearest whole charge
    print protein_charge

    # SET UP LIGAND TOPOLOGY

    print "\nCONSTRUCTING LIGAND TOPOLOGY"

    # Read the ligand into OEMol.
    print "Reading ligand..."
    ligand = readMolecule(ligand_filename, normalize = True)

    # Get formal charge of ligand.
    ligand_charge = formalCharge(ligand)
    print "formal charge is %d" % ligand_charge
    
    # TODO: Choose most appropriate protonation and tautomeric state.
    
    # Write molecule with explicit hydrogens to mol2 file.
    print "Writing ligand mol2 file..."
    ligand_mol2_filename = os.path.abspath(os.path.join(work_path, 'ligand.mol2'))
    writeMolecule(ligand, ligand_mol2_filename)

    # Set substructure name (which will become residue name).
    print "Modifying molecule name..."
    modifySubstructureName(ligand_mol2_filename, 'MOL')

    # Run antechamber to assign GAFF atom types.
    print "Running antechamber..."
    os.chdir(work_path)
    gaff_mol2_filename = os.path.join(work_path, 'ligand.gaff.mol2')
    charge_model = 'bcc'
    command = 'antechamber -i %(ligand_mol2_filename)s -fi mol2 -o ligand.gaff.mol2 -fo mol2 -c %(charge_model)s -nc %(ligand_charge)d > antechamber.out' % vars()
    print command
    output = commands.getoutput(command)    
    os.chdir(current_path)

    # Generate frcmod file for additional GAFF parameters.
    ligand_frcmod_filename = os.path.join(work_path, 'frcmod.ligand')
    output = commands.getoutput('parmchk -i %(gaff_mol2_filename)s -f mol2 -o %(ligand_frcmod_filename)s' % vars())
    print output

    # Run LEaP to generate topology / coordinates.
    ligand_prmtop_filename = os.path.join(work_path,'ligand.prmtop')
    ligand_crd_filename = os.path.join(work_path,'ligand.crd')
    ligand_off_filename = os.path.join(work_path, 'ligand.off')
    
    tleap_input_filename = os.path.join(work_path, 'setup-ligand.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-ligand.leap.out')
    contents = """
# Load GAFF parameters.
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# load ligand
ligand = loadMol2 %(gaff_mol2_filename)s

# check the ligand
check ligand

# report net charge
charge ligand

# save AMBER parameters
saveAmberParm ligand %(ligand_prmtop_filename)s %(ligand_crd_filename)s

# write .off file
saveOff ligand %(ligand_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    ligand_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    ligand_charge = int(round(float(ligand_charge))) # round to nearest whole charge
    print "ligand charge is %d" % ligand_charge    

    # SET UP VACUUM SIMULATION

    # construct pathname for vacuum simulations
    vacuum_path = os.path.join(work_path, 'vacuum')
    if not os.path.exists(vacuum_path):
        os.makedirs(vacuum_path)

    # solvate the ligand
    print "Preparing vacuum ligand with tleap..."
    system_prmtop_filename = os.path.join(vacuum_path,'system.prmtop')
    system_crd_filename = os.path.join(vacuum_path,'system.crd')
    tleap_input_filename = os.path.join(vacuum_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(vacuum_path, 'setup-system.leap.out')
    clearance = 50.0 # clearance in A
    contents = """
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { ligand }
""" % vars()
    # add counterions
    if (ligand_charge != 0):
        nions = abs(ligand_charge)
        if ligand_charge < 0: iontype = 'Na+'
        if ligand_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """

# Create big box.
setBox system centers %(clearance)f

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
    system_prefix = os.path.join(vacuum_path, 'system')
    command = '%(amb2gmx)s --prmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    print command
    output = commands.getoutput(command)
    print output

    # write enlarged box for solvent because LEaP doesn't do it right.
    system_g96_filename = os.path.join(vacuum_path, 'system.g96')
    command = '%s -f %s -o %s -bt octahedron -d %f -center 0.0 0.0 0.0' % (editconf, system_g96_filename, system_g96_filename, clearance)
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
        # add those atoms in the ligand to our list
        if atom['residue'] == 'MOL':
            perturb_atom_indices.append(atom['nr'])
        # add counterions to balance ligand
        if ((atom['residue'] == 'Na+') or (atom['residue'] == 'Cl-')) and (nions < abs(ligand_charge)):
            perturb_atom_indices.append(atom['nr'])
            nions += 1            
    perturbGromacsTopology(system_top_filename, ligand, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = perturb_atom_indices)

    # construct .mdp files
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-vacuum', 'constraints', 'output']))
    mdpfile.write(os.path.join(vacuum_path, 'minimize.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-vacuum', 'constraints', 'thermostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(vacuum_path, 'production.mdp'))
    vacuum_path = os.path.abspath(vacuum_path)

    # write shell script for minimization and production
    contents = """\
#!/bin/tcsh
#BSUB -J %(jobname)s_vacuum
#BSUB -o %(vacuum_path)s/outfile.%%J
#BSUB -e %(vacuum_path)s/errfile.%%J
#BSUB -n 1
#BSUB -M 2000000
#BSUB -W 168:0

source %(GMXRC)s
cd %(vacuum_path)s
setenv RANSEED 1234

# constrained minimize
%(grompp)s -f minimize.mdp -c system.g96 -p system.top -o minimize.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s minimize.tpr -x minimize.xtc -c minimize.g96 -e minimize.edr -g minimize.log

# production
%(grompp)s -f production.mdp -c minimize.g96 -p system.top -o production.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.g96 -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
""" % vars()    
    write_file(os.path.join(vacuum_path, 'run.sh'), contents)

    # PREPARE SOLVATED LIGAND

    print "\nPREPARING SOLVATED LIGAND"

    # create the directory if it doesn't exist
    solvent_path = os.path.join(work_path, 'solvent')
    if not os.path.exists(solvent_path):
        os.makedirs(solvent_path)

    # solvate the ligand
    print "Solvating the ligand with tleap..."
    system_prmtop_filename = os.path.join(solvent_path,'system.prmtop')
    system_crd_filename = os.path.join(solvent_path,'system.crd')
    tleap_input_filename = os.path.join(solvent_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(solvent_path, 'setup-system.leap.out')
    clearance = 10.0 # clearance in A
    contents = """
source leaprc.ff99
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { ligand }
""" % vars()
    # add counterions
    if (ligand_charge != 0):
        nions = abs(ligand_charge)
        if ligand_charge < 0: iontype = 'Na+'
        if ligand_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
solvateOct system TIP3PBOX %(clearance)f

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
    command = '%(amb2gmx)s --prmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    print command
    output = commands.getoutput(command)
    print output

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
        # add those atoms in the ligand to our list
        if atom['residue'] == 'MOL':
            perturb_atom_indices.append(atom['nr'])    
    perturbGromacsTopology(system_top_filename, ligand, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = perturb_atom_indices)

    # set up mdp files
    print "Writing mdp files..."
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'output']))
    mdpfile.write(os.path.join(solvent_path, 'minimize.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'output']))
    mdpfile.setParameter('nsteps', '10000') # 20 ps equilibration
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(solvent_path, 'equilibration.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(solvent_path, 'production.mdp'))

    # write run script
    print "Writing run script..."
    solvent_path = os.path.abspath(solvent_path)
    contents = """\
#!/bin/tcsh
#BSUB -J %(jobname)s_solvent
#BSUB -o %(solvent_path)s/outfile.%%J
#BSUB -e %(solvent_path)s/errfile.%%J
#BSUB -n 1
#BSUB -M 2000000
#BSUB -W 168:0

source %(GMXRC)s
cd %(solvent_path)s
setenv RANSEED 1234

# constrained minimize
%(grompp)s -f minimize.mdp -c system.g96 -p system.top -o minimize.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s minimize.tpr -x minimize.xtc -c minimize.g96 -e minimize.edr -g minimize.log

# equilibration
%(grompp)s -f equilibration.mdp -c minimize.g96 -p system.top -o equilibration.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.g96 -e equilibration.edr -g equilibration.log

# production
%(grompp)s -f production.mdp -c equilibration.g96 -p system.top -o production.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.g96 -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
""" % vars()
    write_file(os.path.join(solvent_path, 'run.sh'), contents)

    # PREPARE SOLVATED COMPLEX

    print "\nPREPARING SOLVATED COMPLEX"

    # create the directory if it doesn't exist
    complex_path = os.path.join(work_path, 'complex')
    if not os.path.exists(complex_path):
        os.makedirs(complex_path)

    # solvate the ligand
    print "Solvating the complex with tleap..."
    system_prmtop_filename = os.path.join(complex_path,'system.prmtop')
    system_crd_filename = os.path.join(complex_path,'system.crd')
    tleap_input_filename = os.path.join(complex_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(complex_path, 'setup-system.leap.out')
    clearance = 10.0 # clearance in A
    contents = """
source leaprc.ff99
source leaprc.gaff
loadamberparams %(parameter_path)s/frcmod_t1p
loadamberparams %(parameter_path)s/frcmod_t2p
loadamberparams %(parameter_path)s/frcmod_y1p
loadamberparams %(parameter_path)s/frcmod_y2p

# load antechamber-generated additional parameters
mods = loadAmberParams %(ligand_frcmod_filename)s

# Load protein.
loadOff %(protein_off_filename)s

# Load ligand.
loadOff %(ligand_off_filename)s

# Create system.
system = combine { protein ligand }
""" % vars()
    # add counterions for ligand
    if (ligand_charge != 0):
        nions = abs(ligand_charge)
        if ligand_charge < 0: iontype = 'Na+'
        if ligand_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions for ligand (will be annihilated with ligand).
addions system %(iontype)s %(nions)d
""" % vars()
    #
    # add counterions for protein
    if (protein_charge != 0):
        nions = abs(protein_charge)
        if protein_charge < 0: iontype = 'Na+'
        if protein_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions for protein.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
solvateOct system TIP3PBOX %(clearance)f

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
    system_prefix = os.path.join(complex_path, 'system')
    command = '%(amb2gmx)s --prmtop %(system_prmtop_filename)s --crd %(system_crd_filename)s --outname %(system_prefix)s' % vars()
    print command
    output = commands.getoutput(command)
    print output

    # set up perturbation
    print "Modifying topology file for perturbation..."
    system_top_filename = os.path.join(complex_path, 'system.top')
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
        # add those atoms in the ligand to our list
        if atom['residue'] == 'MOL':
            perturb_atom_indices.append(atom['nr'])
        # add counterions to balance ligand
        if ((atom['residue'] == 'Na+') or (atom['residue'] == 'Cl-')) and (nions < abs(ligand_charge)):
            perturb_atom_indices.append(atom['nr'])
            nions += 1
            
    perturbGromacsTopology(system_top_filename, ligand, perturb_torsions = True, perturb_vdw = True, perturb_charges = True, perturb_atom_indices = perturb_atom_indices)

    # set up mdp files
    print "Writing mdp files..."
    mdpfile = MdpFile(compose_blocks(['header', 'minimization', 'nonbonded-solvent', 'constraints', 'restraints', 'output']))
    mdpfile.write(os.path.join(complex_path, 'minimize.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'thermostat', 'barostat', 'output']))
    mdpfile.setParameter('nsteps', '10000') # 20 ps equilibration
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(complex_path, 'equilibration.mdp'))

    mdpfile = MdpFile(compose_blocks(['header', 'dynamics', 'nonbonded-solvent', 'constraints', 'restraints', 'thermostat', 'barostat', 'free-energy', 'output']))
    mdpfile.randomizeSeed() # randomize velocities
    mdpfile.write(os.path.join(complex_path, 'production.mdp'))    

    # Copy restraint scripts to run directory
    shutil.copy(compute_angles,complex_path)
    shutil.copy(restraint_topology,complex_path)
    # write run script
    print "Writing run script..."
    complex_path = os.path.abspath(complex_path)
    contents = """\
#!/bin/tcsh
#BSUB -J %(jobname)s_complex
#BSUB -o %(complex_path)s/outfile.%%J
#BSUB -e %(complex_path)s/errfile.%%J
#BSUB -n 1
#BSUB -M 2000000
#BSUB -W 168:0

source %(GMXRC)s
cd %(complex_path)s
setenv RANSEED 1234

# create a tpr file from which to create restraints    
%(grompp)s -f minimize.mdp -c system.g96 -p system.top -o minimize.tpr -maxwarn 10000 -n system.ndx
# Integrate restraints into pre-minimization topfile
python computeangles.py -f system.g96 -s minimize.tpr -d $RANSEED
python restraint_topology.py -n system -p .

# constrained minimize with restraints
%(grompp)s -f minimize.mdp -c system_restr.g96 -p system_restr.top -o minimize.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s minimize.tpr -x minimize.xtc -c minimize.g96 -e minimize.edr -g minimize.log

# equilibration with restraints
%(grompp)s -f equilibration.mdp -c minimize.g96 -p system.top -o equilibration.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s equilibration.tpr -o equilibration.trr -x equilibration.xtc -c equilibration.g96 -e equilibration.edr -g equilibration.log

# create a TPR file, integrate restraints
%(grompp)s -f production.mdp -c equilibration.g96 -p system_restr.top -o production.tpr -maxwarn 10000 -n system.ndx
cp system_restr.top minimize.top
python computeangles.py -f equilibration.g96 -s production.tpr -d $RANSEED
python restraint_topology.py -n equilibration -p .

# production with restraints
%(grompp)s -f production.mdp -c equilibration_restr.g96 -p minimize_restr.top -o production.tpr -maxwarn 10000 -n system.ndx
%(mdrun)s -s production.tpr -o production.trr -x production.xtc -c production.g96 -e production.edr -g production.log

# signal completion
mv run.sh run.sh.done
""" % vars()
    write_file(os.path.join(complex_path, 'run.sh'), contents)

    return

#=============================================================================================
# MAIN
#=============================================================================================

# set up the system
protein_pdb_filename = os.path.join(protein_path, 'j3_1_automodel_mcce_in_mcce_out.pdb')
#protein_pdb_filename = os.path.join(protein_path, 'j3_1_automodel_mcce_out.pdb') # JDC test
ligand_filename = os.path.join(ligand_path, 'jnk.aff-1.sdf')
work_path = 'complex-structures/jnk.aff-1'

# setup system 
setup_system(protein_pdb_filename, ligand_filename, work_path, parameter_path, GMXRC,'jnk_aff_1')


