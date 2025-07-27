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
import time

from util.parser_mode2 import *
from __main__ import *

def bowtie_indexer(TE_fa, transcript_fa, threads):
    ######## Check transcripts fasta IDs if they have repetition
    clock = time()
    print(f"{clock}\tCreating bowtie2 index for TEs...")

    te_fasta = open(TE_fa, 'r')
    i = int(1)

    if check_file(str("index/te_fasta.fa")) == False:
        for line in te_fasta:
            with open(str("index/te_fasta.fa"), 'a') as mod:
                if line.startswith('>'):
                    line = re.sub('\n', '', line)
                    i += 1
                    print(f"{line}_n{i}", file=mod)
                else:
                    line = re.sub('\n', '', line)
                    print(line, file=mod)

    if check_file("index/TE_index.1.bt2") == False:
        subprocess.call(['bowtie2-build', "index/te_fasta.fa", "index/TE_index", '--threads', str(threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(f"TE index has been found!"); print(colored("Skipping...\n", "yellow", attrs=['bold']))



    print(f"{clock}\tCreating bowtie2 index for transcripts...")
    if check_file("index/transcripts_index.1.bt2") == False:
        subprocess.call(['bowtie2-build', transcript_fa, "index/transcripts_index", '--threads', str(threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(f"Transcripts has been index found"); print(colored("Skipping...\n", "yellow", attrs=['bold']))











#
