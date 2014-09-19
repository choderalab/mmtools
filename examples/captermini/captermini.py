#!/usr/bin/env python
#
from mmtools.pipelinetools.pipeline import *

pdbfile = '1YRI.pdb'
seqfile = '1YRI_wtv.seq'
outpdbfile = 'capped.pdb'
thread_model(pdbfile, seqfile, outpdbfile, captermini=True)

