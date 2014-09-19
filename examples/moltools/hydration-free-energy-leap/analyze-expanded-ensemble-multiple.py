#!/usr/bin/python

# Analyze an expanded ensemble simulation conducted with Michael Shirts's modified gromacs 3.1.4 to
# compute free energies.

# NOT FINISHED

#===================================================================================================
# IMPORTS
#===================================================================================================
from numpy import *
from math import *
import MBAR
import timeseries
import commands
import os
import os.path

#===================================================================================================
# CONSTANTS
#===================================================================================================

convert_atmnm3_to_kJmol = 1.01325e5*(1e-09)**3 * 6.02214 * (1e23) # Convert pV from atmospheres*nm^3 into kJ/mol
kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)

#===================================================================================================
# PARAMETERS
#===================================================================================================

temperature = 300.0 # temperature (K)
pressure = 1.0 # pressure (atm)

datafile_basepath = 'molecules/benzamidine/' # directory in which datafiles are stored
datafile_directories  = ['solvent/0'] # directories for phases of the free energy calculation
#datafile_directories  = ['solvent', 'vacuum'] # directories for phases of the free energy calculation
datafile_filename = 'dgdl.xvg' # name of data file in data directories

ndiscard = 40000 # number of snapshots to discard to equilibration TODO take this number from .mdp file (mc-nequil)

#===================================================================================================
# SUBROUTINES
#===================================================================================================

def extract_log_weights(datafile_basepath, datafile_directory):
   """
   """
   
   # Extract weights from log file.
   log_filename = os.path.join(datafile_basepath, datafile_directory, 'production.log')
   print "Reading log weights from %(log_filename)s ..." % vars()
   infile = open(log_filename, 'r')
   log_lines = infile.readlines()
   print "%d lines read." % len(log_lines)
   infile.close()
   g_tk = zeros([T,K], float64)
   t = 0
   while len(log_lines) > 0:
      if log_lines[0].find('MC-lambda information') != -1:
         # eat lines
         log_lines.pop(0)
         log_lines.pop(0)         
         # grab next K lines
         g_k_lines = list()
         for k in range(K):
            line = log_lines.pop(0)
            g_k_lines.append(line)
         # parse lines if we want to keep them
         if (t >= ndiscard):
            for k in range(K):
               line = g_k_lines[k]
               elements = line.split()
               g_tk[t-ndiscard,k] = float(elements[7])
         # increment snapshot counter
         t += 1
         print t
      else:
         # pop top line
         log_lines.pop(0)

   print "g_tk = "
   print g_tk

   return g_tk

#===================================================================================================
# MAIN
#===================================================================================================

# Compute inverse temperature in 1/(kJ/mol)
beta = 1./(kB*temperature)

#===================================================================================================
# Process each phase of the calculation
#===================================================================================================

phases = list() # storage for overall estimate of free energy difference (in kJ/mol)
DeltaF = dict()
dDeltaF = dict()

for datafile_directory in datafile_directories:
   # Construct potential energy datafile name.
   filename = os.path.join(datafile_basepath, datafile_directory, datafile_filename)

   # Read the file.
   infile = open(filename, 'r')
   lines = infile.readlines()
   infile.close()

   # Eliminate header lines.
   while (lines[0][0] == '#'):
      lines.pop(0)

   # Store legend lines
   legend_lines = list()
   while (lines[0][0] == '@'):
      line = lines.pop(0)
      legend_lines.append(line)

   # Discard desired number of snapshots to equilibration.
   for index in range(ndiscard):
      lines.pop(0)

   # DEBUG get rid of last line for now
   lines.pop(len(lines)-1)

   # Determine number of snapshots in file from remaining lines.
   T = len(lines)
   print '%d post-equilibraton snapshots found in dgdl.xvg.' % T

   # Determine number of alchemical intermediates from legend.
   #
   # Legend format looks like this:
   # 
   # @ s0 legend "Current Lambda (fep)"
   # @ s1 legend "dE to Lambda   0.0000(fep)  0.0000(coul)  0.0000(vdw)  0.0000(bonded)  0.0000(restraint)"
   # @ s2 legend "dE to Lambda   0.0000(fep)  0.1000(coul)  0.0000(vdw)  0.0000(bonded)  0.1000(restraint)"
   # ...
   # @ s22 legend "pV Term" [optional]
   K = 0
   for line in legend_lines:
      if line.find('dE to Lambda') != -1:
         K += 1
   print "%d states found." % K
   
   # Get log weights for equilibrium part of calculation.
   log_filename = os.path.join(datafile_basepath, datafile_directory, 'production.log')
   print "Reading log weights from %(log_filename)s ..." % vars()
   g_k_lines = commands.getoutput('grep -A %d "MC-lambda" %s | tail -n %d' % (K+1, log_filename, K)).split('\n')
   g_k = zeros([K], float64)
   # parse lines if we want to keep them
   for k in range(K):
      line = g_k_lines[k]
      elements = line.split()
      g_k[k] = float(elements[7])

   print "g_k = "
   print g_k

   
   # Store current state and reduced potentials at all states.
   print "Extracting reduced potentials..."
   state_t = zeros([T], int32)
   u_tk = zeros([T,K], float64) # u_tk[t,k] is reduced potential of state k from snapshot t
   u_t = zeros([T], float64) # u_t[t] = u_tk[t,k] - g_k[k] where k is current state of system
   for t in range(T):      
      # get state and reduced potentials
      line = lines[t]
      elements = line.split()
      state = int(elements[1]) - 1 # state in range(K)
      current_reduced_potential = float(elements[2])
      for k in range(K):
         u_tk[t,k] = float(elements[3 + k])
      # store
      state_t[t] = state
      u_t[t] = current_reduced_potential - g_k[state]
   print "u_t = "
#   outfile = open('%s.out' % datafile_directory,'w')
#   for t in range(T):
#      outfile.write("%16d %16.8f\n" % (t, u_t[t]))
#   outfile.close()
   g = timeseries.statisticalInefficiency(u_t, u_t)
   print "g = %16.8f" % g
   # compute correlation function
   print "Computing correlation function..."
   C_t = timeseries.normalizedFluctuationCorrelationFunction(u_t, u_t, int(3 * g))
#   outfile = open('corrfun-%s.dat' % datafile_directory, 'w')
#   for t in range(len(C_t)):
#      outfile.write('%8d %16.8f\n' % (t, C_t[t]))
#   outfile.close()
   
   # Test MRS's hypothesis.
   u_t_singlestate = zeros([T], float64)
   for state in range(K):
      # construct timeseries
      Nstate = 0
      for t in range(T):
         if state_t[t] == state:
            #u_t_singlestate[Nstate] = u_tk[t,state]
            u_t_singlestate[Nstate] = u_t[t]
            Nstate += 1

      if Nstate > 0:
         g_state = timeseries.statisticalInefficiency(u_t_singlestate[0:Nstate], u_t_singlestate[0:Nstate])      
         print "state %5d : g = %16.8f, N = %6d" % (state, g_state, Nstate)

   # Analyze timeseries to determine effectively uncorrelated snapshots.
   indices = timeseries.subsampleCorrelatedData(u_t) # indices of uncorrelated samples
   N = len(indices) # number of uncorrelated samples
   print "%d uncorrelated samples of %d snapshots." % (N, T)

   # DEBUG: assume all samples are uncorrelated
#   indices = range(0,T,20)
#   for t in range(T):
#      print "%8d %16.8f" % (t, u_t[t])
#   N = len(indices)
         
   # Count number of uncorrelated samples in each state.
   N_k = zeros(K, int32)
   for n in range(N):
      t = indices[n]
      state = state_t[t]
      N_k[state] += 1
   print "N_k = "
   print N_k
   print "%d uncorrelated samples" % sum(N_k)
   
   # Compute maximum number of samples per state
   N_k_max = max(N_k)

   # Store samples at each state
   u_kln = zeros([K,K,N_k_max], float64) # u_kln[k,l,n] is the reduced potential energy at state l of uncorrelated snapshot n from state k
   N_k = zeros(K, int32) # re-zero N_k
   for n in range(N):
      t = indices[n]
      k = state_t[t]
      u_kln[k,:,N_k[k]] = u_tk[t,:]
      N_k[k] += 1

   # Initialize MBAR (computing free energy estimates, which may take a while)
   print "Computing free energy differences..."
   mbar = MBAR.MBAR(u_kln, N_k, verbose = False, method = 'self-consistent-iteration', CppOptimize = True, maximumIterations = 10000, relative_tolerance = 1.0e-9)
   # mbar = MBAR.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson', CppOptimize = True)   

   # Get matrix of dimensionless free energy differences and uncertainty estimate.
   print "Computing covariance matrix..."
   (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences()
   print Deltaf_ij
   print dDeltaf_ij

   # Print free energy difference
   print "DeltaF = %.3f +- %.3f kcal/mol" % (Deltaf_ij[0,K-1]/beta/4.184, dDeltaf_ij[0,K-1]/beta/4.184)

   # Store
   DeltaF[datafile_directory] = Deltaf_ij[0,K-1]/beta/4.184
   dDeltaF[datafile_directory] = dDeltaf_ij[0,K-1]/beta/4.184   
   
# Compute free energy difference
DeltaF_total = DeltaF['vacuum'] - DeltaF['solvent']
dDeltaF_total = sqrt(dDeltaF['vacuum']**2 + dDeltaF['solvent']**2)

print "total DeltaF = %.3f +- %.3f kcal/mol" % (DeltaF_total, dDeltaF_total)

