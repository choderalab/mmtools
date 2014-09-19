#=============================================================================================
# submsFolders.py
#
# use mmtools to automate setup for the proteins included in our sub-millisecond folders project.
#
# Written 2007-05-22 by
# Gregory R. Bowman <gbowman@stanford.edu>
# Biophysics Program, Stanford University
# Pande group
#=============================================================================================
# TODO:
#=============================================================================================
# Change LOG:
# 5/28/07 - GRB: added code for calling gromacstools
# 6/4/07 - GRB: modified so that each system can have own salt and [salt]
#=============================================================================================
# GLOBAL IMPORTS
import mmtools.mccetools.mcce as mcce
import mmtools.modellertools.modelPDB as modelPDB

# import mmtools.gromacstools.system as system     # the old gromacstools way
from mmtools.gromacstools.System import *

import os
import os.path
import sys
#=============================================================================================

# DEBUG flag: if true then save outputs
DEBUG = True


#----------------------------------
# CLASSES
#----------------------------------

class PipelineProtein:

    def __init__(self, pdb=None, seq=None, salt='NaCl', saltconc=0.100, pH=7.0):
        """Initialize a default data in protein object."""
        
        self.pdbfile  = pdb
        self.seqfile  = seq
        self.salt     = salt         # Molar concentration of salt (ionic strength)
        self.saltconc = saltconc
        self.pH       = pH

        

#----------------------------------
# FUNCTIONS
#----------------------------------
        
 
def shoveItThrough(protein):
    
    # determine base name
    print "PDB File: " + protein.pdbfile
    print "Sequence File: " + protein.seqfile
    extInd = protein.seqfile.rindex(".")
    baseName = protein.seqfile[0:extInd]
    
    # run modelPDB
    if (1):
        modelPDBOut = baseName + "_modelPDB_out.pdb"
        myModel = modelPDB.ModelPDB()
        myModel.makeModel(protein.pdbfile, protein.seqfile, modelPDBOut)
    
    # run mcce
    if (1):
        mcceOut = baseName + "_mcce_out.pdb"
        mcceOut = os.path.abspath(mcceOut)
        prmFile = '../../mccetools/prmfiles/run.prm.quick'
        prmFile = os.path.abspath(prmFile)
        params = mcce.read_paramfile(prmFile)
        mcce.print_prm(params)
        mcce.protonatePDB(modelPDBOut, mcceOut, protein.pH, os.environ['MCCE_LOCATION'], cleanup=True, prmfile=prmFile)
        
    # run gromacs setup
    if (1):
        gromacsOut = baseName + "_final.pdb"
        forcefield = 'ffamber99p'
        # g = system.GromacsSystem(mcceOut, useff=forcefield)    # the old gromacstools way
        g = System(mcceOut, useff=forcefield)
        g.setup.setSaltConditions(protein.salt, protein.saltconc)
        g.setup.set_boxSoluteDistance(0.9)   # <--- ***Greg!!!*** <--- periodic box margin distance, in nanometers 
        thisOutDir = os.path.join(thisdir, baseName)
        print 'Writing equilibration directory to',thisOutDir,'...'
        if os.path.exists(thisOutDir) == False:
            os.mkdir(thisOutDir)
        
        g.prepare(outname=gromacsOut, outdir=thisOutDir, verbose=True, cleanup=False, debug=DEBUG, protocol='racecar2', checkForFatalErrors=True)
    
    # cleanup
    if not DEBUG:
        os.remove(modelPDBOut)
        os.remove(mcceOut)

            

#----------------------------------
# MAIN
#----------------------------------

if __name__ == '__main__':
    

    thisdir = os.getcwd()
    
    proteins = []
    
    # temporarily just test on one pdb
    # files = [["1ENH.pdb", "1ENH_enhd.seq", 'NaCl', 0.150]]
    # files = [["1YRI.pdb", "1YRI_wtv.seq", 'NaCl', 0.150]]
    # files = [["2GB1.pdb", "2GB2_proG.seq", 'NaCl', 0.050]]   # pH 7.5
    files = [["2PTL.pdb", "2PTL_proL.seq", 'NaCl', 0.050]]     # pH 7.0
    # files = [["2F4K.pdb", "2F4K_supv.seq", 'NaCl', 0.050]]     # pH 7.0

    
    #----------------------#
    # PROTEINS             #
    #----------------------#
        
    # Protein L
    proteins.append(  PipelineProtein( pdb="2PTL.pdb", seq="2PTL_proL.seq", salt='NaCl', saltconc=0.050, pH = 7.0) )

    # Protein G
    # proteins.append(  PipelineProtein( pdb="2GB1.pdb", seq="2GB2_proG.seq", salt='NaCl', saltconc=0.050, pH = 7.5) )

    # PSBD
    # proteins.append(  PipelineProtein( pdb="1W4E.pdb", seq="1W4E_psbd.seq", salt='NaCl', saltconc=0.150, pH = 5.5) )

    # BBL
    # proteins.append(  PipelineProtein( pdb="1BAL.pdb", seq="1BAL_bbl.seq", salt='NaCl', saltconc=0.100, pH = 5.3) )

    # NTL9 (wildtype)
    # proteins.append(  PipelineProtein( pdb="1DIV.pdb", seq="1DIV_ntl9.seq", salt='NaCl', saltconc=0.100, pH = 5.4) )

    # NTL9 K12M-G34(dA)
    #proteins.append(  PipelineProtein( pdb="1DIV.pdb", seq="1DIV_ntl9mut.seq", salt='NaCl', saltconc=0.100, pH = 5.4) )

    # Engrailed homeodomain
    #proteins.append(  PipelineProtein( pdb="1ENH.pdb", seq="1ENH_enhd.seq", salt='NaCl', saltconc=0.100, pH = 5.7) )

    # c-myb
    # proteins.append(  PipelineProtein( pdb="1IDY.pdb", seq="1IDY_cmyb.seq", salt='NaCl', saltconc=0.100,  pH = 5.7) )

    # WW Domain 
    # proteins.append(  PipelineProtein( pdb="1I6C.pdb", seq="1I6C_ww.seq", salt='NaCl', saltconc=0.080,  pH = 7.0) )

    # WW Domain mutant
    # proteins.append(  PipelineProtein( pdb="1I6C.pdb", seq="1I6C_wwmut.seq", salt='NaCl', saltconc=0.080,  pH = 7.0) )

    # WW Domain mouse
    # proteins.append(  PipelineProtein( pdb="1E0L.pdb", seq="1E0L_wwmouse.seq", salt='NaCl', saltconc=0.150,  pH = 6.5) )

    # c-myb dP174
    # proteins.append(  PipelineProtein( pdb="1IDY.pdb", seq="1IDY_dP174.seq", salt='NaCl', saltconc=0.100,  pH = 5.7) )

    # villin WT
    # proteins.append(  PipelineProtein( pdb="1YRI.pdb", seq="1YRI_wtv.seq", salt='NaCl', saltconc=0.100,  pH = 7.0) )

    # super-villin 
    # proteins.append(  PipelineProtein( pdb="2F4K.pdb", seq="2F4K_supv.seq", salt='NaCl', saltconc=0.100,  pH = 7.0) )

    # The B-domain of protein A (BDPA)
    # proteins.append(  PipelineProtein( pdb="1BDD.pdb", seq="1BDD_proA.seq", salt='NaCl', saltconc=0.100,  pH = 5.5) )
    

    # THE BUSINESS END #    
    # put each PDB/sequence pair through processing
    for protein in proteins:
        shoveItThrough(protein)
        
