#!/bin/tcsh

# prepare consistent sets of minimization and production files
cat header.mdp minimize.mdp output.mdp nonbonded.mdp free-energy.mdp > minimize-unconstrained.mdp
cat header.mdp minimize.mdp output.mdp nonbonded.mdp constraints.mdp free-energy.mdp > minimize-constrained.mdp
cat header.mdp dynamics.mdp output.mdp nonbonded.mdp bath.mdp constraints.mdp free-energy.mdp > production.mdp

