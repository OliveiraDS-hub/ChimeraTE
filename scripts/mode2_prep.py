import argparse
import datetime
import subprocess
import glob
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
import pybedtools
import contextlib
import sys
import os

clock = time()
print(str(clock) + '\t' + "Creating bowtie2 index for TEs...")

te_fasta = open(str(args.te), 'r')
i = int(1)

for line in te_fasta:
    with open(str(out_dir + '/index/te_fasta.fa'), 'a') as mod:
        if line.startswith('>'):
            line = re.sub('\n', '', line)
            i += 1
            print(line + '_n' + str(i), file=mod)
        else:
            line = re.sub('\n', '', line)
            print(line, file=mod)
    mod.close()
te_fasta.close()

subprocess.call(['bowtie2-build', str(out_dir + '/index/te_fasta.fa'), str(out_dir + '/index/TE_index'), '--threads', str(args.threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
print(colored("Done!", "green", attrs=['bold']))

print(str(clock) + '\t' + "Creating bowtie2 index for transcripts...")
subprocess.call(['bowtie2-build', str(args.transcripts), str(out_dir + '/index/transcripts_index'), '--threads', str(args.threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
print(colored("Done!", "green", attrs=['bold']))










#
