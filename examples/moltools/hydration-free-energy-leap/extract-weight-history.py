#=============================================================================================
# Extract the weight history for all states from production.log file.
#
# Crawls directories 'fecalcs*/*/complex-structures/*/*/production.log'
#
# Written by John D. Chodera <jchodera@gmail.com> 2008-02-13
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

import commands
import os
import os.path
from numpy import *

#=============================================================================================

#=============================================================================================
# PARAMETERS
#=============================================================================================
basepath = "." # base directory to crawl
nlambda = 34 # number of lambda values
final_weights_filename = 'final-weights.out' # master file to write the final weights of each ligand to
avg_weights_filename = 'avg-weights.out' # master file to write the average weights of each ligand to

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
            outfile.write(line + '\n')
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

def format_weights(g_k):
    """Format weights g_k for pretty printing on a line.

    ARGUMENTS
      g_k (numpy array) - list of weights

    RETURNS
      line (string) - the weights printed %10.5 on a line.
    """
    
    line = "%10.5f" % g_k[0]
    for k in range(1,len(g_k)):
        line += " %10.5f" % g_k[k]
    return line

#=============================================================================================
# MAIN
#=============================================================================================

# Make a list of production files
production_file_list = commands.getoutput('ls -1 molecules/*/*/*/production.log').split('\n')

g_k = zeros([nlambda], float64)

final_weights = dict() # final_weights[ligand_name] is a list of numpy vectors containing g_k for each replicate
for production_file in production_file_list:
    print production_file

    # Extract base directory.
    directory = os.path.dirname(production_file)
    print directory
    
    # load contents of file
    production_lines = read_file(production_file)
    nproduction_lines = len(production_lines)
    
    # form g_k.out filename
    weight_filename = os.path.join(directory, 'g_k.out')
    print weight_filename

    # parse weights
    weight_file_contents = ""
    production_line_index = 0
    nsamples = 0
    while (production_line_index < nproduction_lines):
        if production_lines[production_line_index].find('MC-lambda') != -1:
            # advance index to data lines
            production_line_index += 2
            
            # if we don't have a full block, break
            if (production_line_index + nlambda) >= nproduction_lines:
                # incomplete block -- skip it
                break
                        
            # grab the block of MC-lambda            
            block_lines = production_lines[production_line_index:(production_line_index+nlambda)]

            # parse block
            #   N  CoulL   VdwL  BondL  RestL    TempL   Count    F(in kT)  dF(in kT)
            for k in range(nlambda):
                elements = block_lines[k].split()
                lambda_index = int(elements[0])
                log_weight = float(elements[7])
                g_k[k] = log_weight                
                
            # write block to output file lines
            weight_file_contents += format_weights(g_k) + "\n"

            # advance index past block
            production_line_index += nlambda

            # advance sample count
            nsamples += 1
        else:
            # advance index
            production_line_index += 1

    # write weight file
    write_file(weight_filename, weight_file_contents)
    print "%d samples read" % nsamples

    # store final weights
    ligand_name = production_file.split('/')[-3]
    if ligand_name not in final_weights:
        final_weights[ligand_name] = list()
    final_weights[ligand_name].append(g_k.copy()) # append copy of weights to list
        
    # print
    print ""

# write all weights and average weights
final_weights_contents = ""
avg_weights_contents = ""
ligand_name_list = final_weights.keys()
ligand_name_list.sort()
for ligand_name in ligand_name_list:
    g_k_list = final_weights[ligand_name] # list of weights for replicates
    N = len(g_k_list) # number of replicates

    g_k_avg = zeros([nlambda], float64)
    for g_k in g_k_list:
        # accumulate sum
        g_k_avg += g_k
        # write
        final_weights_contents += ligand_name + " " + format_weights(g_k) + "\n"
    # compute average
    g_k_avg /= float(N)
    avg_weights_contents += ligand_name + " " + format_weights(g_k_avg) + "\n"    
    
# write all weights
write_file(final_weights_filename, final_weights_contents)
write_file(avg_weights_filename, avg_weights_contents)


