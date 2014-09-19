#=============================================================================================
# Find the maximal common GAFF substructure among a number of different ligands.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-08
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from openeye.oechem import *
from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
import os
import os.path
from numpy import *

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

#ligand_basepath = '/scratch/users/jchodera/sampl/ligands-parameterized/' # base path for all parameterized ligand files
ligand_basepath = 'ligands-parameterized/' # base path for all parameterized ligand files

# Target ligands
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (8, 10, 11, 15, 22, 26, 32, 38, 52) ] # GOOD - group-1
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (8, 10, 11) ] # GOOD - test
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13, 16, 53) ] # BAD
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13,16,53) ] # BAD
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (13,53) ] # GOOD - group-2
ligand_name_list = [ ('jnk.aff-%d' % index) for index in (49,50) ] # GOOD - group-3
#ligand_name_list = [ ('jnk.aff-%d' % index) for index in (46,60) ] # GOOD - group-4
# intermediate filename (for writing)
intermediate_filename = 'intermediate.mol2'

#=============================================================================================
# SUBROUTINES
#=============================================================================================
#=============================================================================================
# MAIN
#=============================================================================================

# read all ligands
ligands = list()
for ligand_name in ligand_name_list:
    # attempt to load the ligand with GAFF parameters
    ligand = loadGAFFMolecule(ligand_basepath, ligand_name)
    # skip if no ligand found
    if not ligand: continue
    # append to list
    ligands.append(ligand)
print "%d ligands loaded" % len(ligands)

# find common substructure
common_substructure = determineCommonSubstructure(ligands, debug = True)

# find RMS-fit charges
common_substructure = determineMinimumRMSCharges(common_substructure, ligands, debug = True)

# write out molecule
writeMolecule(common_substructure, intermediate_filename)

#ofs = oemolostream()
#ofs.open(intermediate_filename)
#OEWriteConstMolecule(ofs, common_substructure)
#ofs.close()



