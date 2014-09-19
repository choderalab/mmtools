#=============================================================================================
# Example illustrating use of 'moltools' in setting up free energy of annihilation of GAFF phenol in TIP3P solvent.
# 
# Written by John D. Chodera <jchodera@gmail.com> 2008-01-21
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

from mmtools.moltools.ligandtools import *
import commands

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================

#=============================================================================================
# MAIN
#=============================================================================================

# Create a molecule.
# molecule = createMoleculeFromIUPAC('phenol')
# molecule = createMoleculeFromIUPAC('biphenyl')
# molecule = createMoleculeFromIUPAC('4-(2-Methoxyphenoxy) benzenesulfonyl chloride')
molecule = readMolecule('/Users/jchodera/projects/CUP/SAMPL-2008/jnk3/data-from-openeye/jnk.aff/jnk.aff-1.sdf', normalize = True)
   
# Write mol2 file for the molecule.
writeMolecule(molecule, 'phenol.mol2')

# Write GAFF parameters for gromacs, using antechamber to generate AM1-BCC charges.
parameterizeForGromacs(molecule, topology_filename = 'phenol.top', coordinate_filename = 'phenol.gro', charge_model = 'bcc', resname = 'MOL')

# Modify gromacs topology file for alchemical free energy calculation.
perturbGromacsTopology('phenol.top', molecule)

# Convert gromacs topology to an .itp file suitable for inclusion.
top_to_itp('phenol.top', 'phenol.itp', moleculetype = 'phenol')

# Run a free energy calculation.

# enlarge box
command = 'editconf_d -f phenol.gro -o phenol-box.gro -bt octahedron -d 1.0'
print command
output = commands.getoutput(command)
print output

# solvate
topology_contents = """
; amber forcefield
#include "ffamber03.itp"

; my molecule
#include "phenol.itp"

; TIP3P water
#include "ffamber_tip3p.itp"

[ system ]
phenol in water

[ molecules ]
; Compound        nmols
phenol            1
"""
outfile = open('system.top', 'w')
outfile.write(topology_contents)
outfile.close()

command = 'genbox_d -cp phenol-box.gro -cs ffamber_tip3p.gro -o phenol-solvated.gro -p system.top'
print command
output = commands.getoutput(command)
print output

# minmize unconstrained
command = 'grompp_d -f minimization-solvent.mdp -c phenol-solvated.gro -p system.top -o minimization-unconstrained.tpr -maxwarn 10000'
print command
output = commands.getoutput(command)
print output
command = 'mdrun_d -s minimization-unconstrained.tpr -x minimization-unconstrained.xtc -c minimized-unconstrained.gro -e minimization-unconstrained.edr -g minimization-unconstrained.log'
print command
output = commands.getoutput(command)
print output

# minmize constrained
command = 'grompp_d -f minimization-solvent-constrained.mdp -c minimized-unconstrained.gro -p system.top -o minimization-constrained.tpr -maxwarn 10000'
print command
output = commands.getoutput(command)
print output
command = 'mdrun_d -s minimization-constrained.tpr -x minimization-constrained.xtc -c minimized-constrained.gro -e minimization-constrained.edr -g minimization-constrained.log'
print command
output = commands.getoutput(command)
print output

# set up production
command = 'grompp_d -f production-solvent.mdp -c minimized-constrained.gro -p system.top -o production.tpr -maxwarn 10000'
print command
output = commands.getoutput(command)
print output

# run production
command = 'mdrun_d -s production.tpr -o production.trr -x production.xtc -c production.gro -e production.edr -g production.log'
print command
output = commands.getoutput(command)
print output

# make pdb file of trajectory
command = 'echo 0 | trjconv_d -f production.xtc -o production.pdb -s production.tpr'
print command
output = commands.getoutput(command)
print output

