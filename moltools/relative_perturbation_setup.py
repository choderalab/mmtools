#!/bin/env python

#=============================================================================================
# Set up relative alchemical free energy calculations in gromacs using a system prepped in AMBER leap.
#
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-01
# Modified by Michael Zhu <mnz3v@virginia.edu> 2010-09-27
# Rewritten by D. Mobley <dmobley@gmail.com> 2011-03-28 and following
# Main objects of rewrite:
# (1) make it so ligands and protein system (i.e. including ions and possibly ordered waters) come in as leap off files, having been prepared elsewhere.
# (2) Handle input of docking structures of ligands
#=============================================================================================

#CHANGE LOG:
# - D. Mobley, April 2012, placed into svn repository as a starting point for aid in topology/coordinate file setup for relative free energy calculations
# - D. Mobley, May 7, 2012, changed the name as python doesn't like dashes in the file name
# - D. Mobley, May 10, 2012, added exception handling to perturbGromacsTopology to provide a more informative error message in the case where the specified residue name does not exist in the topology file.
# - D. Mobley, May 16, 2012, minor bugfix to a raise of standard error in perturbGromacsTopology
# - D. Mobley, Jan. 7, 2013, added perturbation of mass when atom type is changing.
#                            fixed comments on BAT section when atom types are perturbed but not turned into dummies 
#                            MAJOR BUGFIX: Prior code handled no changes to BAT parameters. These have now been added.                           
# - D. Mobley, Jan. 19, 2013, bugfix -- dummy atom bond and angle torsion parameters for the B state were not being handled correctly (were being perturbed to B state atom types); this is now fixed.
# - D. Mobley, Jan. 22, 2013, bugfix -- transformations wherein one molecule contains an improper torsion which is not present in the other now work (that is, the improper is now turned off in the B state).

#=============================================================================================
# IMPORTS
#=============================================================================================
from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
import pdb
import commands
import os
import os.path
import shutil
import glob
import tempfile
import pipeline
import pickle
import copy



def determine_mutation(ligand1, ligand2, debug = False, atomexpr = OEExprOpts_EqONS, bondexpr = OEExprOpts_BondOrder | OEExprOpts_EqSingleDouble | OEExprOpts_EqAromatic):
    """Take two OEMol molecules corresponding to a ligand and its (precomputed) common substructure with another ligand. Find the best match onto the common substructure and return details for setting up perturbations, such as the retained atoms, their charges, perturbed types, etc.
    
    Written by D. Mobley 1/27/11; motivated by earlier MCSS work by J. Chodera. Updated 4/20/11.
    
    The specific search done here is optimized around bonding patterns; atom matching is extremely loose. This gives great weight to finding rings with matching numbers of atoms, even if the composition of rings is not identical.
    
    Arguments:
     - ligand1: OEMol of molecule 1, corresponding to our ligand
     - ligand2: OEMol of molecule 2, in this case our precalculated common substructure, which DOES NOT NEED FURTHER EDITING
     Optional:
     - debug: False by default; if true, prints extra debugging info
     - atomexpr:  OEExprOpts_EqONS by default, which gives a relatively strict matching when used in combination with the default bondexpr. Option 2 (use with option 2 bond expression) for a much less strict matching:  OEExprOpts_RingMember #Check only ring membership -- atoms match if they share common ring membership 
     - bondexpr:  OEExprOpts_BondOrder | OEExprOpts_EqSingleDouble | OEExprOpts_EqAromatic by default, which gives a relatively strict matching. Option 2 for a much less strict matching: OEExprOpts_ExactBonds | OEExprOpts_EqAromatic | OEExprOpts_EqSingleDouble #CHeck that bonding pattern matches -- either exactly, or counting single and double bonds as the same, or counting aromatic/nonaromatic as the same.

    NOTE on atom and bond expressions: A less strict matching does not necessarily result in more matching atoms, because of conditions on the search.
     Returns:
     - substructure: maximum substructure found
     - numatoms: Number of heavy atoms in match
     - retained_atoms: list of atoms retained in transformation
     - mutated_charges: charges in scaffold that we're mutating to (in dictionary by atom name)
     - mutated_types: atom types in scaffold that we're mutating to (in dictionary by atom name)

    """


    #Initialize common substructure with first ligand
    thismol = ligand1.CreateCopy()

    OEPerceiveChiral( thismol ) #Chirality is not automatically perceived when reading from a file
    OEPerceiveChiral( ligand2 ) #Chirality is not automatically perceived when reading from a file


    initial_substr_size = ligand2.NumAtoms()
    if debug: 
        print "Substructure size %s..." % initial_substr_size


    #Create mcss search from ligand 2

    #Configure search
    #mcss = OEMCSSearch( ligand2, atomexpr, bondexpr )
    mcss = OEMCSSearch( ligand2, atomexpr, bondexpr, OEMCSType_Approximate) #DLM 4/6/11 gives speed increase
    #We still need to do this here as we need to get the matching of our atoms onto the substructure -- but we don't need to EDIT the substructure

    #Set minimum atoms
    mcss.SetMinAtoms( initial_substr_size - 1)
    #Set max matches
    mcss.SetMaxMatches( 5 ) #DLM 4/6/11 modified to 5

    #mcss.SetMCSFunc( OEMCSMaxBondsCompleteCycles() ) #Modify scoring function to prefer completing cycles
    mcss.SetMCSFunc( OEMCSMaxAtomsCompleteCycles() ) #Modify scoring function to prefer completing cycles -- DLM 4/6/11 swapped to this

    if debug: print "Attempting to match %s with substructure..." % (thismol.GetTitle() )

    #Store mutation info (atomname -> charge and atomtype -> typname correspondence for atoms present in substructure)
    mutated_types = dict() #Map of atom types by name
    mutated_charges = dict() # map of charges

    m = OEMol()

    #Perform match and store atom mappings; also track names
    mappings = {}
    name_map = {}
    charge_on_ligand = {}
    ligand_atomtypes = {}
    ct=0
    #map_chiralities = {} #Store chiralities in the map (not original) for reference later
    writeMolecule(thismol, 'test.mol2')
    for match in mcss.Match(thismol, True): #True -> Unique search for common substructure
        ct+=1
        if debug: print "Looking at match number %s..." % ct
        retained_atoms = list() # list of which atoms are not turned into dummies

        #Loop over matched atoms and store atom mappings for later use in checking ring membership
        for matchpair in match.GetAtoms():
            mappings[ matchpair.target ] = matchpair.pattern
            name_map[ matchpair.target.GetName() ] = matchpair.pattern.GetName() #Map names -- scaffold atom names keyed by ligand atom name
            charge = matchpair.pattern.GetPartialCharge() #charge on common substructure
            atomname = matchpair.target.GetName() # name of ligand atom
            charge_on_ligand[ atomname ] = charge
            ligand_atomtypes[ atomname ] = matchpair.pattern.GetType() #Atom type on scaffold
            if debug: print "   Storing charge on ligand atom %s based on scaffold atom %s..." % (atomname, matchpair.pattern.GetName() )

        #Create a match substructure
        m = OEMol()
        OESubsetMol(m, match, True) #Get common substructure
        #This is the main thing we need here

        for atom in  m.GetAtoms():
            thisname = atom.GetName()
            origname = name_map[ thisname ]
            mutated_charges[thisname ] = charge_on_ligand[ thisname ]
            mutated_types[ thisname ]= ligand_atomtypes[ thisname ]
            retained_atoms.append( thisname )

        if debug:
            #Print all retained atoms
            print "Atoms in match and scaffold: "
            atom_index = 1
            
            print "%5s %6s %12s %6s" % ("Index", "Name", "Type", "Name in target molecule")
            for atom in  m.GetAtoms():
                thisname = atom.GetName()
                if name_map.has_key( thisname ):
                    origname = name_map[ thisname ] 
                else: origname = 'none'
                #print "%5d %6s %12s" % (atom_index, atom.GetName(), atom.GetType())
                print "%5d %6s %12s %6s" % (atom_index, thisname, atom.GetType(), origname)
                atom_index += 1
            print ""

        #Only need to consider first match that has all of the atoms in our substructure
        if m.NumAtoms() == initial_substr_size:
            break


    #Ensure the match we are going to return is the correct one
    if m.NumAtoms() == initial_substr_size:
        return m, m.NumAtoms(), retained_atoms, mutated_charges, mutated_types
    else:
        return None, 0, [], {}, {}


#From J. Chodera with modification
def extract_top_section(lines, section):
   """Identify lines associate with a section.

      ARGUMENTS
        lines (list of strings) - the lines in the file
        section (string) - the section name to locate

      RETURNS
        indices (list of integers) - line indices within lines belonging to section
      """

   indices = list()

   nlines = len(lines)
   for start_index in range(nlines):
      # get line
      line = stripcomments(lines[start_index])
      # split into elements
      elements = line.split()
      # see if keyword is matched
      if (len(elements) == 3):
         if (elements[0]=='[') and (elements[1]==section) and (elements[2]==']'):
            # increment counter to start of section data and abort search
            start_index += 1
            break

   # throw an exception if section not found
   if (start_index == nlines):
      raise "Section %(section)s not found." % vars()

   # Locate end of section.
   for end_index in range(start_index, nlines):
      # get line
      line = stripcomments(lines[end_index])
      # split into elements
      elements = line.split()
      # see if keyword is matched
      if (len(elements) == 3):
         if (elements[0]=='['):
           break

   # compute indices of lines in section
   indices = range(start_index, end_index)

   # return these indices
   return indices

def perturbGromacsTopologyToIntermediate(topology_filename, other_ligand_topology, perturbed_resname, ligand, common_substructure, perturb_torsions = True, atomtypes_file = 'ligand_atomtypes.pickle', debug = False):
   """Modify a gromacs topology file to add perturbed-state parameters to transform molecule into a common intermediate.

   ARGUMENTS
     topology_file (string) - the name of the Gromacs topology file to modify
     other_ligand_topology (string) -- name of Gromacs topology containing other ligand. Used in reading atom types to make sure we have a full set of perturbed atom types (i.e. in cases where atom type changes in perturbation)
     perturbed_resname (string) - residue name of residue in topology file to perturb (must be unique)
     ligand (OEMol) - ligand to be perturbed -- must be the same one used to generate the topology file with parameterizeForGromacs(), but with GAFF atom/bond types
     common_substructure (OEMol) - the intermediate to perturb to, with charges assigned and AMBER atom and bondtypes.

   OPTIONAL ARGUMENTS
     perturb_torsions (boolean) - if True, torsions containing at least two atoms in the annihilated region whose central bond is not in an aromatic ring will be turned off in B state (default: True)
     atomtypes_file -- filename of ligand atom types file specifying atomtypes parameters for all atom types. Necessary since MCSS structures occasionally use types which deviate from those in both A and B states (due to bugs in antechamber or oechem) and I need to have a full list of parameters to fill in for those. Default: ligand_atomtypes.pickle
     debug (boolean) -- default False. If True, print extra debug info.

   NOTES
     The common_intermediate MUST have AMBER atom and bondtypes.
     This code currently only handles the special format gromacs topology files produced by amb2gmx.pl -- there are allowed variations in format that are not treated here.
     Note that this code also only handles the first section of each kind found in a gromacs .top file.
     Currently this code handles only dihedrals type 3.
     This code assumes BAT parameters are handled explicitly rather than by atom type.
   
   CHANGELOG:
      MRS: moved to acpype.py
      DLM 3/2011: Moved to much more permissive MCSS matching for consistency with my relative calculation planner. Basically, atoms match if they are in a ring of the same size (regardless of atom/bond type) or if they are in a non-ring structure of the same connectivity map (regardless of atom/bond type and number of connected atoms).
      DLM 4/5/2011: Modified to handle perturbation of atom type when appropriate. 
      DLM 5/20/2011: Modified to make sure perturbed atom types are in atom types section
      DLM 1/7/2013: Various bugfixes, including the MAJOR bugfix that now this handles perturbations to bond, angle, and torsion parameters.

   TODO
     Perhaps this method should be combined with 'parameterizeForGromacs' as an optional second step, ensuring that the correct 'molecule' is used.
     Generalize this code to allow it to operate more generally on a specified moleculetype.

   """

   #Load full dictionary of all potential ligand atom types
   #file = open( atomtypes_file, 'r')
   #ligand_atomtypes = pickle.load( file)
   #file.close()

   #Put all atom types into perturbed atom types -- it can't hurt to have them all there since this is just for the purposes of atom types in the topology file
   if debug: print "Building initial list of perturbed atom types..."
   perturbed_atomtypes = list()
   for atom in ligand.GetAtoms():
       if not atom.GetType() in perturbed_atomtypes:
           perturbed_atomtypes.append( atom.GetType() )
   #Also make sure we have the atom types for the common substructure
   for atom in common_substructure.GetAtoms():
       type = atom.GetType()
       if not type in perturbed_atomtypes:
         perturbed_atomtypes.append( type )

   # find correspondence between common substructure and ligand atoms
   if debug: print "Overlaying molecule onto MCS and determining atom equivalences; this step may take some time..."
   mcss, natoms, retained_atoms, mutated_charges, mutated_types = determine_mutation(ligand, common_substructure, debug = debug)
   if debug:
        print "In resulting match, we found %s matching atoms. Info:" % natoms
        print "   Atoms retained:", retained_atoms
        print "   Mutated types:", mutated_types
 
   if natoms == 0: #Try looser options if no match found
        print "Using looser matching..."
        (mcss, natoms, retained_atoms, mutated_charges, mutated_types) = determine_mutation( ligand, common_substructure, debug = debug, atomexpr = OEExprOpts_RingMember, bondexpr = OEExprOpts_ExactBonds | OEExprOpts_EqAromatic | OEExprOpts_EqSingleDouble )
        if debug:
            print "In resulting match, we found %s matching atoms. Info:" % natoms
            print "   Atoms retained:", retained_atoms
            print "   Mutated types:", mutated_types 

   if natoms==0: #If still zero, stop with warning
        print "Warning: For this match, interpretation of the common substructure was unsuccessful and the perturbation setup will be in error. Pausing..."
        raw_input()

   # Read the contents of the topology file.
   lines = read_file(topology_filename)
   othertop_lines = read_file( other_ligand_topology )


   # Parse [ atomtypes ]
   atomtypes = list() # storage for atom types
   indices = extract_top_section(lines, 'atomtypes')
   existingtypes = list()
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip line if number of elements is less than expected
      if (nelements < 7): continue
      # parse elements
      atomtype = dict()
      atomtype['name'] = elements[0]
      atomtype['bond_type'] = elements[1]
      atomtype['mass'] = float(elements[2])
      atomtype['charge'] = float(elements[3])
      atomtype['ptype'] = elements[4]
      atomtype['sigma'] = float(elements[5])
      atomtype['epsilon'] = float(elements[6])
      existingtypes.append( atomtype['name'])
      # append
      atomtypes.append(atomtype)
      #Double check it is in the perturbed_atomtypes also
      if not atomtype['name'] in perturbed_atomtypes:
        perturbed_atomtypes.append(atomtype['name'])


   insertloc = indices[-1]
   #Parse [atomtypes] from other topology
   indices2 = extract_top_section( othertop_lines, 'atomtypes')
   for index in indices2:
      # extract the line
      line = stripcomments(othertop_lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip line if number of elements is less than expected
      if (nelements < 7): continue
      # parse elements
      atomtype = dict()
      atomtype['name'] = elements[0]
      atomtype['bond_type'] = elements[1]
      atomtype['mass'] = float(elements[2])
      atomtype['charge'] = float(elements[3])
      atomtype['ptype'] = elements[4]
      atomtype['sigma'] = float(elements[5])
      atomtype['epsilon'] = float(elements[6])
      # append
      if not atomtype['name'] in existingtypes:
          atomtypes.append(atomtype)
          existingtypes.append( atomtype['name'] )
          #Form a new line
          line = "%(name)-10s%(bond_type)6s      0.0000  0.0000  A %(sigma)13.5e%(epsilon)13.5e ; copied from B\n" % atomtype
          #Add it after the end of the atom types section we're building
          lines.insert(insertloc, line)



   #Loop over all of the perturbed atom types and make sure that we've got an atom type definition for every type in our topology file now.
   for atype in perturbed_atomtypes:
      if not atype in existingtypes: 
        #If somehow there is an atom type in the substructure that is not in either topology/ligand, then we need to raise an exception and stop
        raise LookupError('Error: When examining atom types present in the perturbation, we found an atom type (%s) for which no atomtypes parameters are specified in either of the ligand topologies. This is an unexpected situation; aborting.' % atype)

   # Augment [ atomtypes ] list with any needed dummy atomtypes with zero LJ well depth.
   indices = extract_top_section(lines, 'atomtypes')
   donetypes = [] #To avoid repeating dummy atoms
   for atomtype in atomtypes:
       if atomtype['name'] in perturbed_atomtypes and not atomtype['name']+'_dummy' in donetypes:
           # make perturbed record
           perturbed_atomtype = atomtype
           # modify name and LJ well depth, leaving other quantities unperturb
           perturbed_atomtype['name'] += '_dummy'
           perturbed_atomtype['epsilon'] = 0.0
           # form the line
           line = "%(name)-10s%(bond_type)6s      0.0000  0.0000  A %(sigma)13.5e%(epsilon)13.5e ; perturbed\n" % perturbed_atomtype
           # insert the new line
           lines.insert(indices[-1], line)
           indices.append(indices[-1]+1)
           #This type is done
           donetypes.append( atomtype['name'] )
    
   #Pre-process [ atoms ] sections of both topology files to build a dictionary of atom numbers to atom types and atom types to masses
   indices = extract_top_section(lines, 'atoms')
   A_atoms_to_types = {}
   B_atoms_to_types = {}
   atomtype_to_mass = {}
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 8): continue
      type = elements[1]
      A_atoms_to_types[ int(elements[0]) ] = type
      mass = float(elements[7])
      if not atomtype_to_mass.has_key(type):
         atomtype_to_mass[type] = mass
   indices = extract_top_section(othertop_lines, 'atoms')
   for index in indices:
      # extract the line
      line = stripcomments(othertop_lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 8): continue
      type = elements[1]
      B_atoms_to_types[ int(elements[0]) ] = type
      mass = float(elements[7])
      if not atomtype_to_mass.has_key(type):
         atomtype_to_mass[type] = mass


   # Process [ atoms ] section
   A_atoms_to_B_types = {} #Storage for storing perturbed type, when applicable (indexed by atom number)
   atoms = list()
   atom_indices = dict()
   perturb_atom_indices = list() # list of atoms that are to be perturbed
   indices = extract_top_section(lines, 'atoms')
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
      atom['massB'] = atom['mass'] #By default no mass perturbation

      # skip if not in speicifed residue
      if not atom['residue'] == perturbed_resname: continue


      # set perturbation type
      atomname = atom['atom']
      # atom type
      if atomname in retained_atoms:
          atom['typeB'] = mutated_types[ atomname ] #Perturb atom type to what it is in the final structure
          atom['chargeB'] = mutated_charges[atomname] # common substructure charge is used
          if atom['type'] <> mutated_types[ atomname]:
            atom['comment'] = 'common substructure; perturbed type'
            atom['massB'] = atomtype_to_mass[ atom['typeB'] ] #Perturb the mass if the type is changing
            #Add this to list of perturbed atoms
            perturb_atom_indices.append(atom['nr'])
          else:
            atom['comment'] = 'common substructure'
      else:
          atom['typeB'] = atom['type'] + "_dummy" # atom is turned into a dummy
          atom['chargeB'] = 0.0 # atom is discharged
          atom['comment'] = 'annihilated'
          perturb_atom_indices.append(atom['nr']) # add atom index to list of atoms that will be perturbed
      #Store perturbed type by atom index
      A_atoms_to_B_types[ atom['nr'] ] = atom['typeB']


      # construct a new line
      line = "%(nr)6d %(type)10s %(resnr)6d %(residue)6s %(atom)6s %(cgnr)6d %(charge)10.5f %(mass)10.6f %(typeB)10s %(chargeB)10.5f %(massB)10.6f ; %(comment)s\n" % atom

      # replace the line
      lines[index] = line

      # store atoms
      atoms.append(atom)

      # store lookup of atom atom names -> atom numbers
      atom_indices[atom['atom']] = atom['nr']

   #Check and see if we actually got any information on atom indices -- if not, there was an error processing the topology file (probably we didn't find the specified molecule) and we should throw an exception
   if len( atom_indices.keys() ) == 0:
      raise StandardError('No valid atoms in specified residue name (%s) found in topology file (%s); probably the wrong residue name was specified. Stopping.' % ( perturbed_resname, topology_filename )) 

   #DLM 1/7/12 editing to add perturbation of bonds, angles, and torsions (aside from deletion of rotatable bonds)
   #Procedure: Look up bond, angle, and torsional parameters by atom type if BAT involves a perturbed atom, otherwise copy A state?
   
   #Pre-process bonds sections to get bond parameters by atom type
   #First process A state
   indices = extract_top_section( lines, 'bonds')
   bond_params = {}
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 5): continue
      types = [ A_atoms_to_types[ int(elements[i])] for i in range(0,2) ]
      key1 = types[0]+'-'+types[1]
      key2 = types[1]+'-'+types[0]
      if not bond_params.has_key( key1) and not bond_params.has_key(key2):
          func = int( elements[2])
          Req = float( elements[3])
          Keq = float( elements[4])
          bond_params[key1] = { 'Req':Req, 'function':func, 'Keq':Keq}         
          bond_params[key2] = bond_params[key1]

   #Then process B state
   indices = extract_top_section( othertop_lines, 'bonds')
   for index in indices:
      # extract the line
      line = stripcomments(othertop_lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 5): continue
      types = [ B_atoms_to_types[ int(elements[i]) ] for i in range(0,2) ]
      key1 = types[0]+'-'+types[1]
      key2 = types[1]+'-'+types[0]
      if not bond_params.has_key( key1) and not bond_params.has_key(key2):
          func = int( elements[2])
          Req = float( elements[3])
          Keq = float( elements[4])
          bond_params[key1] = { 'Req':Req, 'function':func, 'Keq':Keq}         
          bond_params[key2] = bond_params[key1]

   #Pre-process angle sections to get angle parameters by atom type
   #First process A state
   angle_params = {}
   indices = extract_top_section( lines, 'angles')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 6): continue
      types = [ A_atoms_to_types[ int(elements[i]) ] for i in range(0,3) ]
      key1 = types[0]+'-'+types[1]+'-'+types[2]
      key2 = types[2]+'-'+types[1]+'-'+types[0]
      if not angle_params.has_key( key1) and not angle_params.has_key(key2):
          func = int( elements[3])
          theta = float( elements[4])
          cth = float( elements[5])
          angle_params[key1] = { 'cth':cth, 'function':func, 'theta':theta}         
          angle_params[key2] = angle_params[key1]
   #Then process B state
   indices = extract_top_section( othertop_lines, 'angles')
   for index in indices:
      # extract the line
      line = stripcomments(othertop_lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 6): continue
      types = [ B_atoms_to_types[ int(elements[i]) ] for i in range(0,3) ]
      key1 = types[0]+'-'+types[1]+'-'+types[2]
      key2 = types[2]+'-'+types[1]+'-'+types[0]
      if not angle_params.has_key( key1) and not angle_params.has_key(key2):
          func = int( elements[3])
          theta = float( elements[4])
          cth = float( elements[5])
          angle_params[key1] = { 'cth':cth, 'function':func, 'theta':theta}         
          angle_params[key2] = angle_params[key1]
   #Pre-process torsions sections to get torsion parameters by atom type
   #First process A state
   indices = extract_top_section( lines, 'dihedrals')
   dihedral_params = {}
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 11): continue
      types = [ A_atoms_to_types[ int(elements[i]) ] for i in range(0,4) ]
      key1 = types[0]+'-'+types[1]+'-'+types[2]+'-'+types[3]
      key2 = types[3]+'-'+types[1]+'-'+types[2]+'-'+types[0] #If we modify code to allow impropers, this may no longer be true
      if not dihedral_params.has_key( key1) and not dihedral_params.has_key(key2):
          func = int( elements[4])
          if func <> 3: raise "Only [ dihedrals ] function 3 corrently supported."
          #Co-C5 parameters
          C = list()
          for element in elements[5:11]:
            C.append( float(element) )
          dihedral_params[key1] = C
          dihedral_params[key2] = dihedral_params[key1]
   #Then process B state  
   indices = extract_top_section( othertop_lines, 'dihedrals')
   for index in indices:
      # extract the line
      line = stripcomments(othertop_lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 11): continue
      types = [ B_atoms_to_types[ int(elements[i]) ] for i in range(0,4) ]
      key1 = types[0]+'-'+types[1]+'-'+types[2]+'-'+types[3]
      key2 = types[3]+'-'+types[1]+'-'+types[2]+'-'+types[0] #If we modify code to allow impropers, this may no longer be true
      if not dihedral_params.has_key( key1) and not dihedral_params.has_key(key2):
          func = int( elements[4])
          if func <> 3: raise "Only [ dihedrals ] function 3 corrently supported."
          #Co-C5 parameters
          C = list()
          for element in elements[5:11]:
            C.append( float(element) )
          dihedral_params[key1] = C
          dihedral_params[key2] = C


   #NOTE TO SELF: CHECK THAT THIS HANDLES BOTH DIHEDRALS SECTIONS WHEN RELEVANT
   #ALSO THINK ABOUT HOW TO HANDLE DIFFERENT TYPES OF TORSIONS (FIRST PASS: CODE TO HANDLE ONLY TYPE 3)


   # Process [ bonds ] section
   indices = extract_top_section(lines, 'bonds')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 5): continue
      # parse line
      bond = dict()
      bond['i'] = int(elements[0])
      bond['j'] = int(elements[1])
      bond['function'] = int(elements[2])
      bond['Req'] = float(elements[3])
      bond['Keq'] = float(elements[4])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if (bond['i'] not in perturb_atom_indices) and (bond['j'] not in perturb_atom_indices): continue
      #If we get here, it is perturbed, so look up perturbed bond parameters by atom type
      types = [ A_atoms_to_B_types[ int(elements[i]) ] for i in range(0,2) ]
      if 'dummy' in types[0] or 'dummy' in types[1]:
        #If it is being changed into a dummy atom, don't change the parameters
        ReqB = bond['Req']
        KeqB = bond['Keq']
      else:
        key = types[0]+'-'+types[1]
        b = bond_params[key]
        ReqB = b['Req']
        KeqB = b['Keq']

      # construct a new line
      line = "%(i)5d %(j)5d %(function)5d%(Req)12.4e%(Keq)12.4e" % bond
      line += " %(ReqB)12.4e%(KeqB)12.4e ; perturbed\n" % vars() 
      # replace the line
      lines[index] = line

   # Process [ angles ] section
   indices = extract_top_section(lines, 'angles')
   for index in indices:
      # extract the line
      line = stripcomments(lines[index])
      # parse the line
      elements = line.split()
      nelements = len(elements)
      # skip if not all elements found
      if (nelements < 6): continue
      # parse line
      angle = dict()
      angle['i'] = int(elements[0])
      angle['j'] = int(elements[1])
      angle['k'] = int(elements[2])
      angle['function'] = int(elements[3])
      angle['theta'] = float(elements[4])
      angle['cth'] = float(elements[5])
      # skip if an atom range is specified, and this is not among the atoms we're looking for
      if (angle['i'] not in perturb_atom_indices) and (angle['j'] not in perturb_atom_indices) and (angle['k'] not in perturb_atom_indices): continue

      #If we get here, it is perturbed, so look up perturbed angle parameters by atom type
      types = [ A_atoms_to_B_types[ int(elements[i]) ] for i in range(0,3) ]
      #key = types[0].replace('_dummy','')+'-'+types[1].replace('_dummy','')+'-'+types[2].replace('_dummy','')
      key = types[0]+'-'+types[1]+'-'+types[2]
      #If any of the B types isa dummy atom, we want to leave the parameters unperturbed
      if 'dummy' in types[0] or 'dummy' in types[1] or 'dummy' in types[2]:
         thetaB = angle['theta']
         cthB = angle['cth']
      else: #Otherwise look up new B state params
          try: a = angle_params[key]      
          except KeyError:
             #If there is no B state parameter here  
             print line
             print angle_params.keys()
             print key
             raise
          thetaB = a['theta']
          cthB = a['cth']

      # construct a new line
      line = "%(i)5d %(j)5d %(k)5d %(function)5d%(theta)12.4e%(cth)12.4e" % angle
      line += " %(thetaB)12.4e%(cthB)12.4e ; perturbed\n" % vars()
      # replace the line
      lines[index] = line

   # Set rotatable bond torsions in B state to zero, if desired.
   if perturb_torsions:
      # Determine list of rotatable bonds to perturb.
      rotatable_bonds = list()
      for bond in ligand.GetBonds():
         # This test is used because bond.IsRotor() doesn't seem to work correctly (e.g. phenol).
         if (not bond.IsAromatic()) and (bond.GetOrder() == 1) and (bond.GetBgn().GetDegree() > 1) and (bond.GetEnd().GetDegree() > 1):
            i = atom_indices[bond.GetBgn().GetName()]
            j = atom_indices[bond.GetEnd().GetName()]
            if (i < j):
               rotatable_bonds.append( (i,j) )
            else:
               rotatable_bonds.append( (j,i) )

      #print "%d rotatable bond(s) found." % len(rotatable_bonds)

      # Search for [ dihedrals ] section.
      indices = extract_top_section(lines, 'dihedrals') # extract non-blank, non-comment lines
      for index in indices:
         # extract the line
         line = stripcomments(lines[index])
         # parse the line
         elements = line.split()
         nelements = len(elements)
         # skip unrecognized lines
         if (nelements < 11): continue
         # extract dihedral atom indices
         i = int(elements[0])
         j = int(elements[1])
         k = int(elements[2])
         l = int(elements[3])

         # skip if an atom range is specified, and this is not among the atoms we're looking for
         if (i not in perturb_atom_indices) and (j not in perturb_atom_indices) and (k not in perturb_atom_indices) and (l not in perturb_atom_indices): continue

         # function number
         function = int(elements[4])
         if function != 3: raise "Only [dihedrals] function = 3 is supported."
         # C0 - C5 parameters
         C = list()
         for element in elements[5:11]:
            C.append(float(element))

         # reconstruct perturbed line
         line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (i, j, k, l, function, C[0], C[1], C[2], C[3], C[4], C[5])
         if (j,k) in rotatable_bonds:
            # perturb rotatable bonds
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f ; perturbed" % (0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + "\n"
         else:
            # don't perturb 
            line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (C[0], C[1], C[2], C[3], C[4], C[5]) + "\n"

         # replace the line
         lines[index] = line

   #Handle dihedrals section in general, passing over any torsions which are already perturbed
   indices = extract_top_section(lines, 'dihedrals')
   for index in indices:
        #extract line
        line = stripcomments(lines[index])
        elements = line.split()
        nelements = len(elements)
        #Skip unrecognized lines
        if (nelements < 11): continue
        #extract dihedral atom indices
        atom_indices = [ int(elements[i]) for i in range(0,4) ]
        #Obtain perturbed types of these atoms
        types = [ A_atoms_to_B_types[ int(elements[i]) ] for i in range(0,4) ]

        #Figure out if we need to make any changes here; if not, keep going.
        if (atom_indices[0] not in perturb_atom_indices) and (atom_indices[1] not in perturb_atom_indices) and (atom_indices[2] not in perturb_atom_indices) and (atom_indices[3] not in perturb_atom_indices):
                continue

        #Function number
        function = int(elements[4]) 
        if function <> 3: raise "Only [ dihedrals ] function 3 is supported."
        #C0-C5 parameters
        C = list()
        for element in elements[5:11]:
            C.append(float(element))

        #Now, see if any of the atoms are dummies; if they are keep torsions fixed
        if 'dummy' in types[0] or 'dummy' in types[1] or 'dummy' in types[2] or 'dummy' in types[3]:
            types = [ A_atoms_to_types[ int(elements[i])] for i in range(0,4) ] # Use A types
        else:
            types = [ A_atoms_to_B_types[ int(elements[i])] for i in range(0,4) ] # Use B types
        #Otherwise just keep going

        #Look up B state parameters
        try:
           C_B = dihedral_params[ types[0]+'-'+types[1]+'-'+types[2]+'-' + types[3]]
        except KeyError: #If we don't have B state parameters for this -- that is, this torsion occurs in only the molecule we are currently looking at (normally, an improper torsion), then set the B state values to zero. Example: Mutation of ethene to ethane, wherein we are giving up planarity (and an improper) as we do the mutation.
            C_B = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
            
        #ONLY IF NOT ALREADY PERTURBED
        #Reconstruct perturbed line
        if 'perturbed' not in line: #If we're not already perturbing this bond (i.e. turning it to zeros)
           line = "    %-4s %-4s %-4s %-4s %3d%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % (atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3], function, C[0], C[1], C[2], C[3], C[4], C[5])
           # perturb rotatable bonds
           line += " %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f ; perturbed" % (C_B[0], C_B[1], C_B[2], C_B[3], C_B[4], C_B[5]) + "\n"

        # replace the line
        lines[index] = line
   # Replace topology file.
   outfile = open(topology_filename, 'w')
   for line in lines:
      outfile.write(line)
   outfile.close()

   return


def stripcomments(line):
   """Return line with whitespace and comments stripped.
   """
   # strip comments
   index = line.find(';')
   if index > -1:
      line = line[0:index]
   # strip whitespace
   line = line.strip()
   # return stripped line
   return line


