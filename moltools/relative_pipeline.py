#!/bin/env python

#========================================================================================================================
#INFORMATION
#========================================================================================================================
"""This module provides tools relating to setting up relative free energy calculations in GROMACS based on input files in an AMBER format style. Various functionality is provided.

In the initial version, the pipeline here assumes input files are provided as:
- ligand mol2 files in AMBER style (AMBER atom types)
- ligand .prmtop and .crd files (atom names/types matching the mol2 files just mentioned)
- maximal common substructure mol2 file for each planned relative calculation, listing the common substructure for the pair of ligands being considered, also with AMBER atom types
- structure and parameters for the protein the ligand is to be placed into, already in AMBER .prmtop and .crd files.


CHANGELOG:
- Written and adapted from older code early May, 2012, by David Mobley, University of New Orleans
- 5/16/12:
  - Instituted changelog and provided this documentation
  - Migrated input/output stuff to end of file 
  - Wrote setup_ligandpair_relative_topgro, which takes top and gro files for the protein+ligand pair, plus mol2 files for the ligand and common substructure, and does the setup and file bookkeeping
  - Adjusted the above function to automatically read residue names from topology files
  - functionalized more
  - cleaning up intermediate files by placing in a temporary directory which is deleted
- 5/17/12:
  - removed extra debug info relating to mcss search and improved some verbosity
  - added command line interface for a particular use case which comes into play if executed as main
  - improved verbosity


PREREQUISITES:
- GROMACS -- for 'old style' calculations, 4.5.x or 4.0.x should work fine; for new style, 4.6 or similar is recommended
- acpype: http://code.google.com/p/acpype/
- ambertools: http://ambermd.org/AmberTools-get.html (specifically, antechamber should be in your path)
- mmtools from simtk: https://simtk.org/home/mmtools; set your MMTOOLSPATH environment variable. Needs to be mmtools from May 2012 or later.
- OpenEye OEChem and Quacpack (python toolkits), http://www.eyesopen.com
- ffamber ports for gromacs headers/water models: http://ambermd.org/AmberTools-get.html

MAIN FUNCTIONALITY:
- Most applications can potentially be handled using convertPrmCrdToTopGro, ensure_forcefield, and then application of setup_ligandpair_relative_topgro. Other functions present are essentially helper functions utilized by these. 
- If this script is run from the command line, a command line argument parser is provided for one particular use case.

""" 


#========================================================================================================================
#IMPORTS
#========================================================================================================================

from mmtools.moltools.ligandtools import *
import os
import shutil
import commands
import mmtools.moltools.relative_perturbation_setup as pertsetup
import pickle
import tempfile


#========================================================================================================================
#FUNCTION DEFINITIONS
#========================================================================================================================

#Function definitions
def zeroChargesOnDummies( intop, outtop):
    """To break a GROMACS 4.6 perturbation topology into two parts, A doing any charge, LJ, and bonded transformations, and B doing just atomic deletions, this function handles setup of part A from a topology doing the overall calculation. Specifically, it takes an input topology file and locates any atoms being transformed to dummy atoms. The atoms B state type is set to the A state type and any charges are left as zero in the B state. This breaking up is necessary for earlier GROMACS versions where there is only a single lambda value which controls all transformations, rather than GROMACS 4.6 where charge and LJ and other transformations can be controlled separately."""

    #Read input topology
    file = open(intop, 'r')
    text = file.readlines()
    file.close()

    #Extract atoms section and find any dummy atoms, denoted by '_dummy'; switch them to not being dummies in the B state (same type as A state).
    atomlinenrs = extract_section( text, 'atoms' )
    for linenr in atomlinenrs:
        line = text[linenr]
        if '_dummy' in line:
            line = line.replace('_dummy','      ')
            text[linenr] = line
    
    #Write new file
    file = open( outtop, 'w')
    file.writelines(text)
    file.close()


def zeroDummyChargesAndRemovePerturbations( intop, outtop):
    """To break a GROMACS 4.6 perturbation topology into two parts, A doing any charge, LJ, and bonded transformations, and B just doing the LJ transformation for any deletions, this function handles setup of part B from a topology doing the overall calculation. Specificially, it takes an input topology file and locates any dummy atoms, the A state charges for which are zeroed. Then all other perturbations are removed (A state charges set to B state charges; A state types set to B state types except for dummies, and A state bonded parameters set to B state bonded parameters. This breaking up is necessary for earlier GROMACS versions where there is only a single lambda value which controls all transformations, rather than GROMACS 4.6 where charge and LJ and other transformations can be controlled separately.

Handles perturbations to bonds, angles, torsions, and atoms sections. Other sections are ignored."""

    #Read input topology
    file = open( intop, 'r')
    text = file.readlines()
    file.close()

    #Parse atoms section, making appropriate modifications
    atomlinenrs = extract_section( text, 'atoms')
    for linenr in atomlinenrs:
        #Retrieve line
        line = text[linenr]
        linenc, comments = stripcomments(line)
        #Split line
        elements = line.split()        
        nelements = len(elements)
        elements_nc = linenc.split()
        nelements_nc = len( elements_nc )
        #Skip if not all elements found for a perturbed atom
        if nelements_nc < 10: continue

        #Parse line
        atom = dict()
        atom['nr'] = int(elements[0])
        atom['type'] = elements[1] 
        atom['resnr'] = int(elements[2])
        atom['residue'] = elements[3]
        atom['atom'] = elements[4]
        atom['cgnr'] = int(elements[5])                                                                          
        atom['charge'] = float(elements[6])
        atom['mass'] = float(elements[7])
        atom['typeB'] = elements[8]
        atom['chargeB'] = float( elements[9] )
        atom['massB'] = float( elements[10] )
        atom['comment'] = ''
        for elem in elements[12:]:
            atom['comment']+= elem + ' '

        #Make A state charges equal B state charges
        atom['charge'] = atom['chargeB']

        #Make A state type same as B state type unless dummy
        if not '_dummy' in atom['typeB']:
            atom['type'] = atom['typeB']


        #Compose and reinsert line
        line = "%(nr)6d %(type)10s %(resnr)6d %(residue)6s %(atom)6s %(cgnr)6d %(charge)10.5f %(mass)10.6f %(typeB)10s %(chargeB)10.5f %(mass)10.6f ; %(comment)s\n" % atom
        text[linenr] = line
        

    #Next handle bonds section and copy B parameters to A, if present, otherwise keep just A parameters (no change needed)
    indices = extract_section( text, 'bonds')
    for index in indices:
        #Extract the line
        line = text[index]
        (linestripped, comments) = stripcomments(line)
        elements = line.split()
        nelements = len( linestripped.split() ) #Length of line without comments

        #Skip if no B elements found (just keep A parameters)
        if (nelements < 7): continue

        #Parse line
        bond = dict()
        bond['i'] = int( elements[0] )
        bond['j'] = int( elements[1] )
        bond['function'] = int( elements[2] )
        bond['Req'] = float( elements[3] )
        bond[ 'Keq'] = float( elements[4] )
        bond['ReqB'] = float( elements[5] )
        bond[ 'KeqB'] = float( elements[6] )
        bond['comments'] = ''
        for elem in elements[7:]:
            bond['comments'] += elem + ' '
        #Compose a new line
        line = "%(i)5d %(j)5d %(function)5d%(ReqB)12.4e%(KeqB)12.4e %(ReqB)12.4e%(KeqB)12.4e ; %(comments)s\n" % bond
        #Insert line
        text[index] = line

    #Next handle angle section and copy B parameters to A
    indices = extract_section( text, 'angles')
    for index in indices:
        #Extract line with and without comments
        line = text[ index ]
        linenocomments, c = stripcomments( line )

        elements = line.split()
        #Skip line (keep as is) if there is no B state
        if ( len(linenocomments.split() ) < 8): continue

        #Otherwise, get it and make sure A state set to B state
        #Parse line
        angle = dict()
        angle['i'] = int(elements[0] )
        angle['j'] = int(elements[1] )
        angle['k'] = int(elements[2] )
        angle['function'] = int(elements[3] )
        angle['thetaB'] = float(elements[6] )
        angle['cthB'] = float(elements[7] )
        angle['comments'] = ''
        for elem in elements[8:]:
            angle['comments']+= elem + ' '
        #Construct new line
        line = "%(i)5d %(j)5d %(k)5d %(function)5d%(thetaB)12.4e%(cthB)12.4e %(thetaB)12.4e%(cthB)12.4e ; %(comments)s\n" % angle
        #Insert new line
        text[index] = line
    

    #Next handle torsions section and copy B parameters to A
    indices = extract_section( text, 'dihedrals' ) 
    #Note there are multiple possible formats of this section depending on the particular torsional type used

    for index in indices:
        #Extract line with and without comments
        line = text[ index ]
        linenocomments, c = stripcomments( line )
        elements = line.split()
        nelements = len(elements)
        nelements_nocomments = len( linenocomments.split() )

        #Function type is normally element 4; skip if don't have that
        if nelements_nocomments < 4: continue    

        #Else parse
        i = int( elements[0] )
        j = int( elements[1] )
        k = int( elements[2] )
        l = int( elements[3] )
        func = int( elements[4] )

        if func==3:  #Handle type 3 (RBs for AMBER torsions)
            if nelements_nocomments < 16: continue #Skip if not enough data/no perturbation
            C = list()
            for element in elements[5:17]:
                C.append( float (element ))
            comment = ''
            for element in elements[18:]:
                comment+= element+' '

            #Construct new line
            line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f  %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f ; %s \n" % (i, j, k, l, func, C[6], C[7], C[8], C[9], C[10], C[11], C[6], C[7], C[8], C[9], C[10], C[11], comment)

        elif func==9 or func == 1: #Handle type 9, new format, or 1, old format impropers
            if nelements_nocomments < 11: continue #Skip if not enough data/no perturbation
            phase = float( elements[5] )
            kd = float( elements[6] )
            pn = int( elements[7] )
            phaseB = float( elements[8] )
            kdB = float( elements[9] )
            pnB = int( elements[10] )
            comment = ''

            for element in elements[11:]:
                comment+= element + ' '

            #Construct new line
            line = "-4s %-4s %-4s %-4s %3d%12.5f%12.5f%3df%12.5f%12.5f%3d ; %s \n" % (i, j, k, l, func, phaseB, kdB, pnB, phaseB, kdB, pnB, comment ) 



        else:
            raise StandardError('[ dihedrals ] function type %s not supported...' % function ) 

        #Insert line
        text[index] = line



    #Write out resulting topology
    file = open(outtop, 'w')
    file.writelines( text)
    file.close()


def rungmx( cmd, checkForFatalErrors = True, verbose = False):
    """Execute a specified gromacs command on the command line, storing the output and checking it for fatal errors. Returns the output as well. Raises an exception when it encounters an error.

Optional arguments:
- maxwarn: Default, 50 -- adjust maximum warnings tolerated before GROMACS gives a fatal error."""

    if verbose: print '>>',cmd
    output=commands.getoutput(cmd)
    if verbose: print output

    if checkForFatalErrors:
        if output.count('Fatal error') > 0:
            print "There was a fatal error executing the following command: %s, run in %s" % (cmd, os.getcwd())
            raise RuntimeError('Error running (in directory %s) the gromacs command %s...' % (os.getcwd(),cmd))

    return output


def getForcefields():
    """Obtain and return list of force fields with codes to allow user to select desired force field. Returns forcefields, forcefieldcodes. Requires gmxlib environment variable to be set."""

    gmxlib = os.getenv('GMXLIB') 
    if gmxlib == None:
        raise RuntimeError('Error: Please set GMXLIB environment variable to point to location of force fields before running.')
    else:
        try:
            fin = open( os.path.join( gmxlib, 'FF.dat'), 'r')
        except IOError:
            raise RuntimeError('Error: FF.dat cannot be found in GMXLIB=%s...' % gmxlib)
            sys.exit(1)
        forcefields = {}
        forcefieldCodes = {}
        lines = fin.readlines()
        numff = int(lines.pop(0))
        for i in range(0, numff):
            fields = lines[i].split()
            key = fields.pop(0)
            forcefields[ key ] = string.joinfields(fields)
            forcefieldCodes[ key ] = str(i)
        fin.close()
        return forcefields, forcefieldCodes


def writeBasicMinimizationMDP( filename, nsteps = 10 ):
    """Write a really basic minimization mdp file for steepest descents minimization out to specified filename. Performs specified number of steps of minimization."""

    input="""
; RUN CONTROL PARAMETERS =
integrator               = steep
nsteps                   =%(nsteps)s
; Output frequency for energies to log file and energy file =
nstlog                   = 1
nstenergy                = 1
; ENERGY MINIMIZATION OPTIONS =
; Force tolerance and initial step-size =
emtol                    = 100
emstep                   = 0.01
; Max number of iterations in relax_shells =
niter                    = 20
; Number of correction steps to use for L-BFGS minimization
nbfgscorr                = 10
; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency =
nstlist                  = 1
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.0
domain-decomposition     = no
; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = pme
rcoulomb                 = 1.0
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure =
DispCorr                  = AllEnerPres
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; OPTIONS FOR BONDS     =
constraints              = none

free_energy = no
""" % vars()

    file = open( filename, 'w')
    file.write( input )
    file.close()


def getSolventGroupFromGenion( tprfile ):
    """Get number of the 'SOL' group from genion and return it. Requires tpr file as input as genion has that requirement."""

    startdir = os.getcwd()
    tempdir = tempfile.mkdtemp() #Temporary directory
    shutil.copy( tprfile, os.path.join( tempdir, 'tmp.tpr') ) 
    os.chdir(tempdir)
    genion_out = commands.getoutput('echo q | genion -s tmp.tpr -np 1')
    genion_lines = genion_out.split('\n')
    groupstr = None
    for line in genion_lines:
        if line[0:5] == 'Group':
            fields = line.replace(')',' ').replace('(',' ').split()
            if fields[2].strip() == 'SOL':
                groupstr = fields[1]
                break
    if groupstr == None:
        print genion_lines
        raise RuntimeError('Error in getting solvent group from genion on tpr file %s (working directory %s); no solvent group found.' % (tprfile, os.getcwd()) )
   

    shutil.rmtree(tempdir)
    os.chdir(startdir)
    return int(groupstr)

def getNetCharge( top, gro, mdpfile, maxwarn = None):
    """Take a top and gro file and an mdp file for an run of user choice (which will not be performed) and run grompp, which calculates the net charge of the system. Extract the net charge and return it. Optional maxwarn argument can be set to an integer; if it is, this is passed to gromacs as the maxwarn option."""

    startdir = os.getcwd()
    #Copy to a temporary directory to avoid generating output files in the working directory
    tempdir = tempfile.mkdtemp()
    shutil.copy( top, os.path.join( tempdir, 'tmp.top') )
    shutil.copy( gro, os.path.join( tempdir, 'tmp.gro') )
    shutil.copy( mdpfile, os.path.join( tempdir, 'tmp.mdp') )
    os.chdir(tempdir)


    #Run grompp and extract output
    if maxwarn:
        rungmxcmd = 'grompp -f tmp.mdp -p tmp.top -c tmp.gro -maxwarn %(maxwarn)s' % vars()
    else:
        rungmxcmd = 'grompp -f tmp.mdp -p tmp.top -c tmp.gro'
    gromppout = commands.getoutput( rungmxcmd )
    
    #Find total charge in output
    totalcharge = None
    for line in gromppout.split('\n'):
        if 'total charge' in line and not 'State B' in line:
            line = line.replace('System has non-zero total charge:','')
            charge = float( line.split()[0].strip() )
            totalcharge = int( round(charge) )
            break #Break at first occurrence, to avoid going on to, for example, state B charge if present.
    if not totalcharge:
        #If the net charge is already zero, it will not be printed. We can attempt to read it from the topology file, as this will often be tracked there. The following code will do so, but assumes the topology file contains only a single [ atoms ] block which is followed by a [ bonds ] block; it reads the last line of the section preceding the [ bonds ] block and atempts to find 'qtot' at the end of the line. 
        file = open('tmp.top', 'r')
        text = file.readlines()
        file.close()
        bondsstart = text.index('[ bonds ]\n')
        lastatomline = bondsstart-1
        while len(text[ lastatomline]) < 20 or text[lastatomline].find( ';' ) < 2: #While the line is too short or starts with a comment, keep reading back
            lastatomline -=1
        #Now extract qtot from line if present
        if 'qtot' in text[ lastatomline ]:
            lastsec = text[ lastatomline ].split('qtot')[1]
            qtot = lastsec.split()[0].strip()
            totalcharge = float(qtot)
        else: 
            #Attempt to get it by summing the charges in the topology file's atoms section
            indices = extract_section( text, 'atoms' )
            totalcharge = None
            for index in indices:
                line, comments = stripcomments( text[index] )
                elements = line.split()
                if len( elements ) >  7:
                    if totalcharge == None: totalcharge = 0.0
                    totalcharge += float(elements[6])
            if totalcharge == None:
                #Otherwise give up and we'll raise an exception below
                print "No total charge found in topology file."

        if totalcharge == None: #If we still can't find the total charge
            print "ERROR while calculating total charge; Work done in %s..." % os.getcwd()
            raise RuntimeError('Error: Could not read total charge from grompp output: %s' % gromppout )

    #Clean up
    os.chdir( startdir ) 
    shutil.rmtree(tempdir)

    #Return
    return totalcharge

def do_editconf( ingro, outgro, size=1.2, box='dodecahedron', vector = None, verbose = False):
    """Use GROMACS editconf to edit the simulation box for the input gro file and put the new box in the specified output gro file. 
    Arguments:
        - ingro: Input gro file
        - outgro: Output gro file
    Option arguments:
        - size: In nm, specifying distance from solute (centered) to nearest box edge. Default 1.2 nm
        - box: A valid gromacs box specificier, such as 'dodecahedron' (default)
        - vector: A gromacs box vector (alternative to size); must be three-element list. Overrides size and passes this to gromacs.
        - verbose: Set verbosity. Default false
    Returns:
        - output: Output of editconf
    Note exception will be raised if editconf encounters fatal error.
    """

    if not vector:
        gmxcmd = 'editconf -f %(ingro)s -o %(outgro)s -bt %(box)s -d %(size)s' % vars() 
    else:
        gmxcmd = 'editconf -f %s -o %s -bt %s -box %s %s %s' % (ingro, outgro, box, vector[0], vector[1], vector[2])
    output = rungmx( gmxcmd, verbose = verbose)

    return output

def make_big_vdwradii( targetpath ):
    """Make big vdwradii.dat in targetpath with modified (big, 0.3 nm) vdW radii to be used when solvating a protein to make sure no waters get buried in the protein where they shouldn't be, as is sometimes the case when using gromacs tools for solvating."""

    file = open( os.path.join( targetpath, 'vdwradii.dat'), 'w')
    text="""; Very approximate VanderWaals radii
; only used for drawing atoms as balls or for calculating atomic overlap.
; longest matches are used
; '???' or '*' matches any residue name
; 'AAA' matches any protein residue name
; MODIFIED TO USE BIG VDW RADII TO PREVENT WATERS BEING PUT IN THE PROTEIN WHERE WE DON'T WANT THEM. DLM
???  C     0.3
???  F     0.3
???  H     0.3
???  N     0.3
???  O     0.3
???  P     0.3
???  S     0.3
???  LP1   0
???  LP2   0
SOL  H     0.04
SOL  O     0.105
WAT  H     0.04
WAT  O     0.105
GLY  MN1   0
GLY  MN2   0
ALA  MCB1  0
ALA  MCB2  0                                                                            
VAL  MCG1  0                                                                            
VAL  MCG2  0                                                                            
ILE  MCG1  0                                                                            
ILE  MCG2  0                                                                            
ILE  MCD1  0                                                                            
ILE  MCD2  0                                                                            
LEU  MCD1  0                                                                            
LEU  MCD2  0                                                                            
MET  MCE1  0                                                                            
MET  MCE2  0                                                                            
TRP  MTRP1 0                                                                            
TRP  MTRP2 0
THR  MCG1  0
THR  MCG2  0
LYSH MNZ1  0                                                                            
LYSH MNZ2  0                                                                            
""" 
    file.writelines(text)
    file.close()


def do_solvate( ingro, outgro, inouttop, water = 'ffamber_tip3p.gro', verbose = False, use_big_vdwradii = False, FF = 'ffamber99sb'):
    """Use genbox to add water to simulation box, except for regions within vdwradii.dat of protein/ligand.
    Arguments:
    - ingro: Input gro file
    - outgro: Output gro file
    - inouttop: Input/output topology -- edited only by adding solvent line and ensuring water model matches that specified in water argument
    Optional:
    - Water: Water coordinate file (forming basis for water model to be included in top); default ffamber_tip3p.gro
    - use_big_vdwradii (boolean): Use modified larger vdw radii designed to prevent solvation water from getting stuck in tiny buried cavities in protein where it may have trouble reaching equilibrium
    - FF: Which force field to use in solvated topologies, if not already present (default: ffamber99sb)
    Returns:
    - output: Output of genbox.
    - nwaters: Number of waters added. 

    Notes: Attempts to edit topology file to ensure you are using the correct water model, based in the water model name specified as an argument. This edit will fail if the water itp file has a different name from the gro file used for solvating. If the topology file you are running this on does not include either spc.itp or your water model already, the water model include will be added at the top of the file. 
    Also ensures that the 

    """

    #Copy top/gro files to a temporary working directory where we will do this setup before copying back. This is necessary because of issues with 3.3.x genbox, requiring input/output topology to be in working directory
    tempdir = tempfile.mkdtemp()
    shutil.copy(inouttop, os.path.join(tempdir, 'tmp.top') )
    shutil.copy(ingro, os.path.join( tempdir, 'tmp.gro') )
    startdir = os.getcwd()
    os.chdir( tempdir )

    if use_big_vdwradii: 
       make_big_vdwradii('.') 

    #Ensure forcefield
    ensure_forcefield( 'tmp.top', 'tmp.top', FF = FF )

    #Run genbox 
    if verbose: print "Running in temporary directory %s..." % tempdir
    cmd = 'genbox -cp tmp.gro -p tmp.top -o out.gro -cs %(water)s' % vars()
    output = rungmx( cmd, verbose = verbose)
    if verbose: print output

    #Compile name of itp file and edit topology to include this
    wateritp = water.replace('.gro','.itp')
    #Edit topology file to replace spc.itp with this
    file = open( 'tmp.top', 'r')
    text = file.readlines()
    file.close()
    try:
        idx = text.index('#include "spc.itp"\n')
        text[idx ] = '#include "%(wateritp)s"\n' % vars()
    except ValueError:
        #Check if it already includes correct water model
        try:
           #If this works then we leave as is
           idx = text.index('#include "%(wateritp)s"\n' % vars() )
        except ValueError:
            #If it is still not found, add the includejust before the SYSTEM section
            idx = 0
            while '[ system ]' not in text[idx]: #Find system line
                idx+=1
            text[idx] = '#include "%(wateritp)s"\n' % vars() + text[idx]
    file = open( 'tmp.top', 'w')
    file.writelines( text )
    file.close()
        
    #Parse output and figure out the number of waters, typically SOL
    tmp = output.split('\n')
    tmp.reverse()
    for line in tmp:
        if 'Number of SOL ' in line: break
    if not 'Number of SOL' in line: nwaters = None #Not found
    else:
        nwaters = int( line.split()[-1] )

    #Copy files back and cleanup
    os.chdir(startdir)
    shutil.copy( os.path.join(tempdir, 'tmp.top'), inouttop )
    shutil.copy( os.path.join(tempdir, 'out.gro'), outgro )
    shutil.rmtree( tempdir)

    return output, nwaters

def add_counterions( intop, ingro, outtop, outgro, negion = 'Cl-', posion ='Na+', negionchg = -1, posionchg = 1, verbose = False, maxwarn = None ):
    """Add counterions to make the system net neutral.
    Generates a temporary mdp file to run genion; reads grompp output to determine net charge of system and adds positive ornegative ions to make system net neutral.
    ARGUMENTS:
        - intop: Name of input topology
        - ingro: Input coordinate file
        - outtop: Output topology containing ions
        - ougro: Output gro file containing ions
    OPTIONAL:
        - negion: Specifier for negative ion type. Default: Cl-
        - posion: Specifier for positive ion type. Default: Na+
        - negionchg: Net charge of negative ion (default -1)
        - posiionchg: Charge of positive ion (default 1)
        - maxwarn: Specify maximum number of warnings for grompp (default: use grompp default). Use this option if there is a warning that you want to ignore...
    To do:
        - Add optional argument to allow system to be solvated at specified ionic strength, in addition to being net neutral.

    Returns:
        - numpos, numneg: number of positive and negative ions actually added."""

    startdir = os.getcwd()
    #Make a temporary directory to work in
    tempdir = tempfile.mkdtemp()

    #Make sure we have absolute paths for output top and gro
    outtop = os.path.abspath( outtop )
    outgro = os.path.abspath( outgro )

    #Copy files to temp dir
    shutil.copy(intop, os.path.join( tempdir, 'in.top' ) )
    shutil.copy(ingro, os.path.join( tempdir, 'in.gro' ) )

    #Change directory there
    os.chdir( tempdir ) 

    #Write basic mdp file
    writeBasicMinimizationMDP( 'tmp.mdp' )

    #Generate tpr
    if maxwarn:
        rungmxcmd = 'grompp -f tmp.mdp -p in.top -c in.gro -o tmp.tpr -maxwarn %(maxwarn)s' % vars()
    else:
        rungmxcmd = 'grompp -f tmp.mdp -p in.top -c in.gro -o tmp.tpr'
    output = rungmx( rungmxcmd )
    if verbose: print output

    #Determine number of solvent group
    solgroup = getSolventGroupFromGenion( 'tmp.tpr' )

    #Determine net charge of system (for now we are just adding enough counterions to neutralize)
    netcharge = getNetCharge( 'in.top', 'in.gro', 'tmp.mdp' )
    if verbose: print "Net charge is %.2g" % netcharge

    #Calculate number and type of ions to add
    tolerance = 1e-2 #Tolerance for zeroing the system net charge
    if netcharge==0 or abs(netcharge) - tolerance < 0:
        numion = 0
        numpos = 0
        numneg = 0
    elif netcharge > tolerance:
        #Check to see if we can neutralize the system
        if abs(netcharge % abs(negionchg)) > tolerance:
            raise RuntimeError("Error: Cannot neutralize the system using ions of charge %s because the netcharge of %d is not divisible by that number; fraction is %.2g." % ( negionchg, netcharge, abs(netcharge/negionchg)))
        else:
            numion = int( floor( netcharge/float(abs(negionchg)) ) )
            numneg = numion
            numpos = 0
    else:
        #Check to see if we can neutralize the system
        if abs(netcharge % abs(posionchg)) > tolerance:
            print netcharge % posionchg
            raise RuntimeError("Error: Cannot neutralize the system using ions of charge %s because the netcharge of %d is not divisible by that number; fraction is %.2g." % ( negionchg, netcharge, abs(netcharge/negionchg)))
        else:
            numion = int( floor( abs(netcharge)/float(posionchg) ) )
            numpos = numion
            numneg = 0

    #Now we are prepared to add (numion) counterions of type (addion) using genion, replacing (solgroup) 
    if numion > 0:
        if verbose: print "Running genion to add %s negative and %s positive ions." % (numneg, numpos)
        genioncmd = 'echo %(solgroup)s | genion -s tmp.tpr -o out.gro -p in.top -np %(numpos)s -nn %(numneg)s -nname %(negion)s -pname %(posion)s' % vars()
        output = rungmx( genioncmd, verbose = verbose )
        if verbose: print output
        os.chdir(startdir)
        shutil.copy( os.path.join( tempdir, 'in.top'), outtop )
        shutil.copy( os.path.join( tempdir, 'out.gro'), outgro )
    else:
        if verbose: print "System is already neutral; no ions added."
        shutil.copy( os.path.join( tempdir, 'in.top'), outtop )
        shutil.copy( os.path.join( tempdir, 'in.gro'), outgro )

    if verbose: print "Finished adding ions and final topology/coordinates created."


    #Parse resulting topology file and make sure it includes ions.itp; if not, add it prior to the system section
    file = open( outtop, 'r')
    text = file.readlines()
    file.close()
    found = False
    for line in text:
        if 'ions.itp' in line:
            found = True
            break
    if not found:
        idx = text.index('[ system ]\n')
        text[idx] = '#include "ions.itp"\n\n' + text[idx] 
        file = open( outtop, 'w')
        file.writelines(text)
        file.close()

    #Clean up
    shutil.rmtree( tempdir )
    os.chdir(startdir)

    return numpos, numneg

def ensure_forcefield( intop, outtop, FF = 'ffamber99sb'):
    """Open a topology file, and check to ensure that includes the desired forcefield itp file. If not, remove any [ defaults ] section (which will be provided by the FF) and include the forcefield itp. Useful when working with files set up by acpypi -- these need to have a water model included in order to work, and most water models require a force field included in order for them to work.

ARGUMENTS:
   - intop: Input topology
   - outtop: Output topology
OPTIONAL:
   - FF: STring corresponding to desired force field; default ffamber99sb.

Limitations:
  - If you use this on a topology file that already includes a DIFFERENT forcefield, the result will be a topolgoy file including two forcefields.
"""

    file = open (intop, 'r')
    text= file.readlines()
    file.close()


    #Check if force field is included somewhere
    found = False
    for line in text:
        if FF in line:
            found = True
    #If not, add it after any comments at the top
    if not found:
        idx = 0
        while text[idx].find(';')==0:
            idx+=1
        text[idx] = '#include "%s.itp"\n\n' % FF + text[idx]

    #Remove any defaults section
    found = False
    foundidx = None
    endidx = None
    for (idx, line) in enumerate(text):
        if '[ defaults ]' in line:
            foundidx = idx
            found = True
        #If we've already found the defaults section, find location of start of next section
        #Assumes next section can be found by looking for a bracket at the start of a line 
        elif found and '[' in line:
            #Only store if we didn't already store
            if endidx == None:
                endidx = idx 
    #Now remove defaults section
    
    if found:
        text = text[0:foundidx] + text[endidx:]

        
    #Write results
    file = open( outtop, 'w')
    file.writelines(text)
    file.close()


def ensure_solvent( intop, outtop, solventinclude = 'ffamber_tip3p.itp'):
    """Reads a topology file and checks to ensure it includes the appropriate solvent include file; if not, this is inserted before the system section. Also can be used to check for any other include files which should come right before the [ system ] section, such as ions.itp."""
    file = open(intop, 'r')
    text = file.readlines()
    file.close()

    found = False
    for line in text:
        if solventinclude in line:
            found = True
            break
    if not found:
        for idx, line in enumerate(text):
            if '[ system ]' in stripcomments(line):
                text[idx] = '#include "%s"\n\n' % solventinclude + text[idx]

                break
    file = open( outtop, 'w')
    file.writelines(text)
    file.close()


def parse_input_topgro_names( name ):
    """Helper function. Take a file name/pair of file names as used by setup_ligandpair_relative_topgro and check for their existence with assertions, sorting out and returning the actual filenames. Used when we are unsure whether a variable contains either (a) a single prefix specifying the prefix of a .top and .gro file, or (b) a list/tuple containing a pair of top/gro names."""

    #Check whether we're working with a list or prefix
    if not os.path.isfile(name[0]):
        #If the first entry is not a name, then it is probably a prefix
        names = (name + '.top', name + '.gro')
        for n in names:
            assert os.path.isfile(n), "No such input file %s..." % n
        return names 
    else:
        names = name
        for n in names:
            assert os.path.isfile(n), "No such input file %s..." % n

        return names


def get_ligand_resname_from_topology( topfile ):
    """Parse specified topology file which is expected to contain only a single residue name. Return residue name as specified in the atoms section."""

    file = open( topfile, 'r')
    text = file.readlines()
    file.close()
    indices = extract_section(text, 'atoms')
    for linenr in indices:
        line, comments = stripcomments(text[linenr])
        elements = line.split()
        if len(elements) > 5:
            resname = elements[3]
            return resname
 
         

def setup_ligandpair_relative_topgro( input_protein_topgro, input_ligand_topgro, input_ligand_mol2, input_mcss_mol2, output_complex_topgro, output_ligand_topgro, debug = False, verbose = True, oldstyle = False, add_after_resnum = None):

    """Takes a pair of ligands and a protein as GROMACS .top and .gro files, and additionally .mol2 files for the ligands and maximal common substructure, and sets up GROMACS .top and .gro files for transformation of each ligand into the common substructure, both in solution and in gas phase.  

ARGUMENTS:
(Note: All arguments can be provided as paths or just as file names, in which case they are assumed to be in the current directory.)
- input_protein_topgro: Prefix of protein top/gro files as a single entry if they share the same name (i.e. 'protein' for 'protein.top' and 'protein.gro'), or a list/tuple specifying the full filename of first the topology file and then the gro file
- input_ligand_topgro: A list for the pair of ligands being considered, where the first entry is for the first ligand and the second entry is for the second ligand. Each entry may be either a file prefix (prior to .top and .gro) or a tuple/list specifying the full filename of the .top and .gro
- input_ligand_mol2: A list for the pair of ligands being considered, specifying the full paths 
- input_mcss_mol2: filename (including path if desired) of mcss file for common substructure of pair of ligands being considered
- output_ligand_topgro: Top/gro file prefixes for the complex and ligands being considered, where the first entry is the output name for the first ligand and the second entry is the output name for the second ligand. These names are used for the solvated complex.
- output_ligand_topgro: Top/gro file prefixes for the output of the ligands being considered, where the first entry is the output name for the first ligand and the second entry is the output name for the second ligand. These names are used for the solvated ligands.

OPTIONAL:
- Debug, boolean specifying whether to print debug info. Default false.
- Verbose, boolean specifying whether to print progress info and other verbosity. Default True.
- oldstyle: Boolean as to whether to generate topologies appropriate for GROMACS 4.5.x and earlier (additional intermediate topologies). If these are to be generated, their names will be created by adding '_oldstyle_LJ' and '_oldstyle_charge' before the extension of the filenames specified in output_complex_topgro and output_ligand_topgro.
- add_after_resnum: Specify at what point in coordinate file for protein ligand is to be added, by residue number. (Residue number it's to be added after). Default: None, meaning at end of coordinates. Useful to add ligand before any ordered waters, for example.

RETURNS:

EXAMPLES:
- setup_ligandpair_relative_topgro( '2vci', ['30d', '40f'], ['30d_amber.mol2', '40f_amber.mol2'], 'mcss.mol2' ['complex_30d_to_40f_pert_solvated', 'complex_40f_to_30d_pert_solvated'], ['30d_to_40f_pert_solvated', '40f_to_30d_pert_solvated'] ) 
- setup_ligandpair_relative_topgro( '2vci', [('30d.top', '30d.gro'), ('40f.top', '40f.gro')], ['../setup/30d_amber.mol2', '../setup/40f_amber.mol2'], '../setupmcss.mol2' [ ('output/complex_30d_to_40f_pert_solvated.top', 'output/complex_30d_to_40f_pert_solvated.gro)', ('output/complex_40f_to_30d_pert_solvated.gro', 'output/complex_40f_to_30d_pert_solvated.gro')], ['30d_to_40f_pert_solvated', '40f_to_30d_pert_solvated'] ) 
"""

    #Do bookkeeping on file names to make sure we have things the way we need
    #For protein input, check whether we're working with a list or a prefix and get names
    input_protein_top, input_protein_gro = parse_input_topgro_names( input_protein_topgro)
    #For ligands, check whether we're working with list or prefix and get names   
    pair_tops = [] 
    pair_gros = []
    for entry in input_ligand_topgro:
        t, g = parse_input_topgro_names( entry)
        pair_tops.append(t)
        pair_gros.append(g) 
    
    #For output, get full names of output tops and gros by starting with input names and checking for .gro and .top
    output_complex_tops = []
    output_ligand_tops = []
    output_complex_gros = []
    output_ligand_gros = []
    for entry in output_complex_topgro:
        if '.top' in entry[0]:
            output_complex_tops.append( entry[0] )
            output_complex_gros.append( entry[1] )
        else:
            output_complex_tops.append( entry + '.top')
            output_complex_gros.append( entry + '.gro' )
    for entry in output_ligand_topgro:
        if '.top' in entry[0]:
            output_ligand_tops.append( entry[0] )
            output_ligand_gros.append( entry[0] )
        else:
            output_ligand_tops.append( entry + '.top')
            output_ligand_gros.append( entry + '.gro')

    #Build a list of residue names for ligands by parsing their topology files
    resnames = []
    for filename in pair_tops:
        resnames.append( get_ligand_resname_from_topology( filename ) )


    #FILE BOOKKEEPING:
    #Load ligands
    if verbose: print "   Reading ligands and mcss..."
    ligandMols = [ readMolecule(file) for file in input_ligand_mol2 ]
    #Load mcss
    mcss_mol = readMolecule( input_mcss_mol2 )

    #Generate temporary directory for intermediate files
    tempdir = tempfile.mkdtemp()
    if verbose: print "   Generating intermediate files in temporary directory %s..." % tempdir


    #GENERATE PERTURBED TOPOLOGIES
    #This tool requires an OEMol of the intermediate to perturb to, with AMBER atom and bond types (?).
    #Apply to generate topologies for mutation of both ligands
    for idx, top in enumerate(pair_tops):
        #FIgure out name of other ligand topology
        if idx==0: othertop = pair_tops[1]
        else: othertop = pair_tops[0]

        #Copy ligand top/gro to temporary intermediates for editing
        intermediate_topology = os.path.join( tempdir, 'ligand%s.top' % idx)
        shutil.copy( top, intermediate_topology )
        intermediate_gro = os.path.join( tempdir, 'ligand%s.gro' % idx)
        shutil.copy( pair_gros[idx], intermediate_gro )

        #Perturb topology
        if verbose: print "   Generating perturbed topologies; this step may take some time due to the MCSS search."
        pertsetup.perturbGromacsTopologyToIntermediate( intermediate_topology, othertop, resnames[idx] , ligandMols[idx], mcss_mol, debug = debug) 

        #ONLY FOR GROMACS 4.5.x and earlier -- break the perturbation into two stages, first turning off the charges on any deleted atoms, then secondarily making a transformation to dummies for any deleted atoms. Any nondeleted atoms will have their charges/types modified in the first step, and any changes to bonded parameters will also be made in that step.
        if oldstyle: #If we want to generate old style topologies...
            #Generate the new file names
            intermediate_OC_topology = intermediate_topology.replace('.top', '_oldstyle_charge.top')
            intermediate_LJ_topology = intermediate_topology.replace('.top', '_oldstyle_LJ.top')

            if verbose: print "   Generating old style GROMACS 4.5.x and earlier topologies for ligand %s/2..." % (idx+1)
            #Generate new topologies for step A, charge/bonded transformation
            zeroChargesOnDummies( intermediate_topology, intermediate_OC_topology) 

            #Generate new topology for step B, LJ deletion
            zeroDummyChargesAndRemovePerturbations( intermediate_topology, intermediate_LJ_topology)

        #Generate top/gro files for the complex by adding the ligand to the protein. Then solvate and add counterions 
        if  verbose: print "   Final system prep -- generating complex, solvating, adding counterions for ligand %s/2..." % (idx+1)
        
        #Generate combined gro file by adding ligand after protein
        if verbose: print "      Generating complex for ligand %s/2..." % (idx+1)

        #Generate intermediate filenames for the complex
        intermediate_complex_topology = os.path.join( tempdir, 'complex%s.top' % idx )
        intermediate_complex_gro = os.path.join( tempdir, 'complex%s.gro' % idx )

        #Generate combined gro file
        if not add_after_resnum:
            add_ligand_to_gro( input_protein_gro, pair_gros[idx], intermediate_complex_gro, resname = resnames[idx] )
        else:
            add_ligand_to_gro( input_protein_gro, pair_gros[idx], intermediate_complex_gro, resname = resnames[idx], add_after_resnum = add_after_resnum)
        #Generate combined topology file for perturbation by adding ligand to topology
        add_ligand_to_topology( input_protein_top, intermediate_topology, intermediate_complex_topology)

        #SOLVATE
        #First generate box
        #For complex
        if verbose: print "      Solvating..."
        #Solvate complex in dodecahedral simulation box with distance at least 1.2 nm from the protein+ligand to the nearest box edge
        #Distance to box edge can perhaps be adjusted downwards depending on target cutoffs and size of protein
        do_editconf( intermediate_complex_gro, intermediate_complex_gro, size = 1.2)
        #For ligand in solution
        do_editconf( intermediate_gro, intermediate_gro, size = 1.2)
    
        #Get ready to solvate
        #First generate intermediate topologies
        intermediate_solvated_topology = os.path.join( tempdir, 'ligand%s_solvated.top' % idx)
        intermediate_complex_solvated_topology = os.path.join( tempdir, 'complex%s_solvated.top' % idx )
        shutil.copy( intermediate_topology, intermediate_solvated_topology )
        shutil.copy( intermediate_complex_topology, intermediate_complex_solvated_topology )

        #Solvate the complex, using big vdw radii for atoms in the system to ensure waters aren't placed close to overlapping/at relatively tiny buried sites in the protein
        output, nwaters_complex = do_solvate( intermediate_complex_gro, intermediate_complex_gro, intermediate_complex_solvated_topology, use_big_vdwradii = True, verbose = verbose ) 

 
        #Now for ligand in solution
        #Solvate, using larger vdW radii for atoms in the system to ensure waters aren't placed close to overlapping
        output, nwaters_ligand = do_solvate( intermediate_gro, intermediate_gro,  intermediate_solvated_topology, use_big_vdwradii = True, verbose = verbose ) 
        
        #ADD COUNTERIONS -- often this applies only to the complex, but here we'll do both in the off chance that the ligand carries a net charge (note that if that is the case, one may need to worry about other complications as well).
        #Handle the complex
        if verbose: print "      Adding counterions..."
        numpos_complex, numneg_complex = add_counterions( intermediate_complex_solvated_topology, intermediate_complex_gro, intermediate_complex_solvated_topology, intermediate_complex_gro, maxwarn = 50, verbose = verbose)
        #Handle the ligand
        numpos_ligand, numneg_ligand = add_counterions( intermediate_solvated_topology, intermediate_gro, intermediate_solvated_topology, intermediate_gro, maxwarn = 50, verbose = verbose)

        #If we are generating old style topologies as well, we need to do the topology setup for the intermediate calculations also
        if oldstyle:
            if verbose: print "      Generating topologies and solvating for old style setup as well..."
            for topname in ['_oldstyle_charge', '_oldstyle_LJ']:
                #Generate intermediate files for complex
                intermediate_complex_topology_oldstyle = os.path.join( tempdir, 'complex%s%s.top' % (idx, topname))

                #Generate topologies for the complex
                add_ligand_to_topology( input_protein_top, os.path.join( tempdir, 'ligand%s%s.top' % (idx, topname)), intermediate_complex_topology_oldstyle )

                #Edit topologies -- both for the complex and for the ligand in solution -- to add the appropriate number of waters and counterions; not necessary to do this for gro files as we can use the same gro files
                #First handle the complex
                file = open( intermediate_complex_topology_oldstyle, 'a')
                #Compute number of solvent molecules to add
                nwaters = nwaters_complex
                if numpos_complex <> None:
                    nwaters-= numpos_complex
                if numneg_complex <> None:
                    nwaters-= numneg_complex
                #Add waters
                file.write('SOL \t %d\n' % nwaters )
                #Determine whether to add positive ions and if so, add them (assuming sodium)
                if numpos_complex > 0 :
                    file.write('NA+ \t %s\n' % numpos_complex)
                #Determine whether to add negative ions and if so, add them (assuming chloride)
                if numneg_complex > 0:         
                    file.write('Cl- \t %s\n' % numneg_complex)
                file.close()
                #Make sure FF for solvent is included
                ensure_solvent( intermediate_complex_topology_oldstyle, intermediate_complex_topology_oldstyle) #By default uses ffamber_tip3p
                #Now also ensure ions.itp is included
                ensure_solvent(  intermediate_complex_topology_oldstyle, intermediate_complex_topology_oldstyle, solventinclude = 'ions.itp' )

                #Next handle the ligand
                tmptop = os.path.join( tempdir, 'ligand%s%s.top' % (idx, topname) )
                file = open( tmptop, 'a')
                #Compute number of solvent molecules to add
                nwaters = nwaters_ligand
                if numpos_ligand <> None:
                    nwaters -= numpos_ligand
                if numneg_ligand <> None:
                    nwaters -= numneg_ligand
                #Add waters
                file.write('SOL \t %d\n' % nwaters )
                #Determine whether to add positive ions and if so, add them (assuming sodium)
                if numpos_ligand > 0 :
                    file.write('NA+ \t %s\n' % numpos_ligand)
                #Determine whether to add negative ions and if so, add them (assuming chloride)
                if numneg_ligand > 0:         
                    file.write('Cl- \t %s\n' % numneg_ligand)
                file.close() 
                #Make sure FF for solvent is included
                ensure_solvent( tmptop, tmptop )

                #Copy any old style files back as appropriate
                #First generate final file names
                couttop = output_complex_tops[idx].replace('.top', '%s.top' % topname )
                coutgro = output_complex_gros[idx].replace('.gro', '%s.gro' % topname )
                louttop = output_ligand_tops[idx].replace('.top', '%s.top' % topname )
                loutgro = output_ligand_gros[idx].replace('.gro', '%s.gro' % topname )
                #Copy files
                shutil.copy( intermediate_complex_topology_oldstyle, couttop )
                shutil.copy( intermediate_complex_gro, coutgro )
                shutil.copy( tmptop, louttop )
                shutil.copy( intermediate_gro, loutgro )

        #Now, copy back files for the general case 
        shutil.copy( intermediate_complex_solvated_topology, output_complex_tops[idx] )
        shutil.copy( intermediate_complex_gro, output_complex_gros[idx] )
        shutil.copy( intermediate_solvated_topology, output_ligand_tops[idx] )
        shutil.copy( intermediate_gro, output_ligand_gros[idx] )

    #Clean up temporary directory
    shutil.rmtree( tempdir )




#========================================================================================================================
#INPUT/OUTPUT STUFF
#========================================================================================================================

#If executing as main, take input to drive...
if __name__=='__main__':

    #Use command line input via option parser
    from optparse import OptionParser
    usage = """usage: %prog -l (ligfile) -p (protname)
where (ligfile) is the name of a tab/space delimited text file, where the pair of ligands to be considered are the first two entries on each line, and the filename for their maximal common substructure (in mol2 format) is the third entry on each line. These prefixes are used to find (ligname).crd and (ligname).top, parameter and coordinate files, and (ligname)_amber.mol2, the .mol2 structure files. (protname) is the prefix of the name of the protein (AMBER .top and .crd files) which the ligand will be combined with for binding calculations.

PURPOSE:
Takes pairs of ligands for planned relative binding free energy calculations and plans transformations to mutate each ligand of each pair to the (provided) maximal common substructure. The starting point is prepared amber parameter and coordinate files for each ligand and the protein, as well as .mol2 (assumed, but .sdf would likely also work) files for each ligand and the maximal common substructure for each pair. All files must use AMBER-style atom types and the atom names and types in the provided .mol2 files must match those in the amber parameter files.

REQUIRED ARGUMENTS:
- '-l' or '--ligandfile' -- specifies filename of ligand file described above
- '-p' or '--protein' -- specifies protein prefix as described above
OPTIONAL ARGUMENTS:
- '-i' or '--inputfile_directory' -- specifies location searched for protname.top and .crd, the ligand input .top and .crd files, and the ligand input .mol2 file, as well as the mcss file. 
- '-o' or '--oldstyle': Also set up topoologies for transformations in GROMACS 4.5.x and earlier; output names will have 'oldstyle' in them. 
- '-q' or '--quiet': Turn off verbosity. 

EXAMPLE:
python relative_pipeline.py -l ligfile.txt -p 2vci ../setup --oldstyle

"""
    parser = OptionParser(usage)
    parser.add_option('-l', '--ligandfile', help = 'Text file (tab/space delimited) specifying the ligand pairs to be considered and their corresponding mcss files, where the first column is the name of ligand1, the second column is the name of ligand2, and the third column si the filename of the mcss file, assumed to be in inputfile_directory. Ligand names are the prefixes used for the top/gro/mol2 files (mol2 files will be found under (name)_amber.mol2.', dest = 'ligfile', metavar = 'FILE')
    parser.add_option('-p', '--protein', help = "Base name for protein the ligands are to be placed into, i.e. for 2vci.crd and .top (prmtop) specify '2vci' here. Assumed to be in inputfile_directory", dest = 'protein_basename', metavar = 'PREFIX')
    parser.add_option('-i', '--inputfile_directory', help = 'Specify inputfile_directory, the directory searched for input ligand files, mcss file, and protein input. Default: Current working directory.', default = os.getcwd(), dest = 'inputfile_directory', metavar = 'DIR')
    parser.add_option('-o', '--oldstyle', help = 'Flag specifying whether or not to also output topologies suitable for free energy calculations using GROMACS 4.5.x and earlier (oldstyle) rather than just for new style GROMACS (4.6 and later). Default: False. Setting the flag turns this on.', dest = 'oldstyle', action = 'store_true', default = False )
    parser.add_option('-q', '--quiet', help = 'Turn off verbosity (be more quiet).', dest = 'verbose', action = 'store_false', default = True )
    #Parse args
    (options, args ) = parser.parse_args()

    #Check for required args
    if not options.ligfile or not os.path.isfile( options.ligfile ):
        parser.error('Error: Please enter valid name of ligand file specifying pairs of ligands to be considered and their mcss files.')
    if not options.protein_basename or not os.path.isfile( os.path.join( options.inputfile_directory, options.protein_basename+'.top')) or not os.path.isfile( os.path.join( options.inputfile_directory, options.protein_basename +'.crd' ) ):
        parser.error('Error: You are required to specify a valid protein basename that points to valid .top and .crd files in the specified input file directory. Files not found.' )

    #Bring some things into current scope
    protein_basename = options.protein_basename
    inputfile_directory = options.inputfile_directory

    #Parse the ligand file to get the required info
    file = open( options.ligfile, 'r')
    text = file.readlines()
    file.close()
    ligandlist = [] #To be a list of all ligands we are working with; code below will look for mol2 files 'name'+'_amber.mol2' where name is the names specified here. Also will look for parameter and coordinate files 'name'.top and 'name'.crd in AMBER format.
    ligandpairs = [] #List containing pairs of ligands we want to perturb between; each entry should be a pair
    mcss = [] #List of MCSS structure files (mol2 format, sdf probably also will work)
    for line in text:
        elements = line.split()
        if len(elements) >= 3 and not line[0]=='#':
            ligandlist.append( elements[0] )
            ligandlist.append( elements[1] )
            ligandpairs.append( [ elements[0], elements[1] ] )
            mcss.append( elements[2] )

    #Set up old style (gromacs 4.5.x) perturbations as well? Boolean
    oldstyle = options.oldstyle

    
    #Convert to GROMACS format and ensure forcefield is present
    for mol in ligandlist + [ protein_basename ]:
        #Convert top and gro -- use GROMACS 4.0.x and 4.5.x style topologies (rather than converting torsions in a way that will only work with 4.5.x and later)
        convertPrmCrdToTopGro( os.path.join( inputfile_directory, mol + '.top'), os.path.join( inputfile_directory, mol +'.crd' ), mol + '.top', mol+'.gro', style = 'Old' )
        #Ensure forcefield is present in topology -- this doesn't affect the parameters (which are all explicitly specified) but will make sure standard atom types like those used in water models are defined
        ensure_forcefield( mol + '.top', mol + '.top', FF = 'ffamber99sb' )


    #Loop over ligand pairs and actually run the topology setup
    input_protein_topgro = os.path.join( protein_basename )
    for idx, pair in enumerate(ligandpairs):

        #First generate desired input/output names
        input_ligand_topgro = [  nm  for nm in pair ]
        input_ligand_mol2 = [ os.path.join( inputfile_directory, nm+'_amber.mol2') for nm in pair ]
        output_complex_topgro = [ 'complex_%s_to_%s_pert' % (pair[0], pair[1]), 'complex_%s_to_%s_pert' % (pair[1], pair[0] ) ] 
        output_ligand_topgro = [ 'ligand_%s_to_%s_pert' % (pair[0], pair[1]), 'ligand_%s_to_%s_pert' % (pair[1], pair[0] ) ]

        #Apply function to do setup 
        setup_ligandpair_relative_topgro( input_protein_topgro, input_ligand_topgro, input_ligand_mol2, os.path.join( inputfile_directory, mcss[idx]), output_complex_topgro, output_ligand_topgro, oldstyle = True, verbose = options.verbose ) 
         
