#!/usr/bin/env python

#=============================================================================================
# Example of resolvating a GROMACS structure with the same number of water molecules
# as a template *.gro file.
#
# In this example resolvate a random-coil conformation of NTL9(1-39) with the same number of 
# waters from a template *.gro of the solvated native structure
# 
# Written by Vincent Voelz <vvoelz@stanford.edu> 2008-04-21
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

import sys, os, commands, glob, tempfile

# GLOBAL IMPORTS
from mmtools.pipelinetools import *
from mmtools.pdbtools import *

import mmtools.modellertools.modelPDB as modelPDB
import mmtools.mccetools.mcce as mcce

from mmtools.gromacstools.MdpFile import *
from mmtools.gromacstools.System import *
from mmtools.gromacstools.MdpFile import *


# GLOBAL VARIABLES

VERBOSE = True


#=============================================================================================
# molecular system setup
#=============================================================================================

pdbfile = '1DIV.pdb'             # The native PDB structure of the N-Terminal domain of ribosomal protein L9 (NTL9)
seqfile = '1DIV_ntl9-39.seq'     # The seqeuence of NTL9, res 1-39.

coilpdb = 'coilconf39mer.pdb'    # A random-coil polyglycine PDB that we wish to resolvate
groTemplate = '1DIV_ntl9-39.gro' # The solvated native structure of NTL9(1-39), which will be used as
                                 # as a template to resolvate the random-coil sturcture 
				 
outdir = 'ntl9_coil_resolvated'
			
# create a PipelineProtein object
protein = pipeline.PipelineProtein( pdb=pdbfile, seq=seqfile, salt='NaCl', saltconc=0.015,  pH = 4.8)
if VERBOSE:
    print
    print '=== created PipelineProtein ==='
    print protein.info()
    
    
# Fixed settings
forcefield = 'ffamber03'
protocol = 'racecar3'
cap = True
cutoff = 1.0    
box = 'small'

# create an output directory for the pipeline work
if os.path.exists(outdir) == False:
    os.mkdir(outdir)
protein.setup(outdir)

# Read in the sequence from the protein.seqfile
fseq = open(protein.seqfile,'r')
sequence = fseq.readline().strip()
fseq.close()

# read in the maximum box dimension from the groTemplate
g = GromacsStructureFromGrofile(groTemplate)
boxstr = g.getboxstring()
maxboxdim = 0.0
for field in boxstr.split():
    if (float(field) > maxboxdim):
	     maxboxdim = float(field)
    print 'boxstr', boxstr
    print 'maxboxdim', maxboxdim


# Thread sequence onto this model PDB
print pipeline.thread_model(coilpdb, protein.seqfile, protein.modelPDBout, captermini=cap)

# convert the pdb into a *.gro file
g = GromacsStructureFromPDB(protein.modelPDBout)

# set the (correctly-protonated) sequence of the template groTemplate
gtemplate = GromacsStructureFromGrofile(groTemplate)
gtemplate_seq = gtemplate.atoms.sequence()
print 'gtemplate.atoms.sequence', gtemplate_seq

# map the correctly-protonated template sequence onto the coil groStructure
g = mapSequenceOntoGromacsStructure(gtemplate_seq, g)

# copy the box dimensions from the template to the coil groStructure
g.boxstring = boxstr

# write the coil.gro (no solvent yet) to a tmpdir
tmpdir = tempfile.mkdtemp()   # temporary directory for unfolded pdbs 
coilgrofile = os.path.join(tmpdir,'coil')
print 'Writing', coilgrofile, '...'
g.write( coilgrofile )   # gromacstools.GromacsData forces us to save this using only the basename
coilgrofile = coilgrofile + '.gro'     # correct to reflect the full name

# write a single molecule of water *.gro to the tmpdir
growater = os.path.join(tmpdir,'ffamber_tip3p_single.gro')
fout = open( growater, 'w')
fout.write("""TIP3P 
3
    1SOL     OW    1   1.336   1.583   0.888  0.1910 -0.3183 -0.5328
    1SOL    HW1    2   1.399   1.553   0.953 -0.9817 -0.4598  0.5822
    1SOL    HW2    3   1.287   1.652   0.933 -2.0169 -1.8033 -0.5086
   2.48706   2.48706   2.48706

""")
fout.close()

# prepare a resolvated grofile from the tempolate grofile -- this will correct the protonation state too
s = System(coilgrofile, finalOutputDir=outdir, finalOutputName=protein.basename, useff=forcefield)   # a gromacstools.System.System() object
s.prepareResolvatedFromTemplate( coilgrofile, groTemplate, growater, protocol=protocol)



