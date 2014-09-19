#=============================================================================================
# Example illustrating use of 'moltools' in computing hydration free energies of small molecules.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-01-23
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.moltools.ligandtools import *
import commands
import os    

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

ligand_path = 'ligands' # path with input ligands
output_path = 'gromacs' # path to put gromace parameter

#=============================================================================================
# MAIN
#=============================================================================================

# Make a list of ligand filenames
ligand_filenames = os.listdir(ligand_path)

# parameterize each ligand
for filename in ligand_filenames:
    # construct full pathname
    fullpath = os.path.join(ligand_path, filename)
    print "Processing %(fullpath)s ..." % vars()

    # read molecule
    molecule = readMolecule(fullpath, normalize = True)    

    # get name by truncating filename
    (prefix, suffix) = os.path.splitext(filename)

    # DEBUG
    formal_charge = formalCharge(molecule)    
    print "%64s %8d" % (prefix, formal_charge)
    
    # Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
#    parameterizeForGromacs(molecule, topology_filename = os.path.join(output_path, prefix + '.top'), coordinate_filename = os.path.join(output_path, prefix + '.gro'), charge_model = 'bcc', resname = 'LIG')
    
    # Convert .top file to .itp file.
#    top_to_itp(topfile = os.path.join(output_path, prefix + '.top'), outputitp = os.path.join(output_path, prefix + '.itp'), moleculetype = 'ligand')
    
