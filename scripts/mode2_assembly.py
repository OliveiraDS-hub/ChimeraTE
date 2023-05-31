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
import shutil

def transcriptome_assembly():
    from __main__ import aln_dir
    from __main__ import trinity_out
    from __main__ import mate1
    from __main__ import mate2

    if str(args.strand) == "rf-stranded":
        strand = "RF"
    else:
        strand = "FR"

    clock = time()
    print(str(clock) + '\t' + "Performing transcriptome assembly... It may take a while")
    subprocess.call(['Trinity', '--seqType', 'fq', '--SS_lib_type', str(strand), '--max_memory', str(str(args.ram) + 'G'), '--left', str(mate1), '--right', str(mate2), '--CPU', str(args.threads), '--output', str(trinity_out)], stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
    print(colored("Done!", "green", attrs=['bold']))
    shutil.rmtree(str(trinity_out + '/chrysalis'), ignore_errors=False, onerror=None)
    shutil.rmtree(str(trinity_out + '/insilico_read_normalization'), ignore_errors=False, onerror=None)
    shutil.rmtree(str(trinity_out + '/read_partitions'), ignore_errors=False, onerror=None)
    os.remove(str(trinity_out + '/scaffolding_entries.sam'))
    os.remove(str(trinity_out + '/jellyfish.kmers.25.asm.fa'))
    os.remove(str(trinity_out + '/both.fa'))

    clock = time()
    print(str(clock) + '\t' + "Creating bowtie2 index with assembled transcripts...")
    subprocess.call(['bowtie2-build', str(trinity_out + '/Trinity.fasta'), str(trinity_out + '/trinity_assembly'), '--threads', str(args.threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(str(clock) + '\t' + "Alignment with assembled transcripts...")
    subprocess.call(['bowtie2', '-x', str(trinity_out + '/trinity_assembly'), '-1', str(mate1), '-2', str(mate2), '-p', str(args.threads), '-S', str(trinity_out + '/transcripts_aln.sam')], stdout=subprocess.DEVNULL)
    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(str(clock) + '\t' + "Converting alignment to bed...")
    with open(str(trinity_out + '/transcripts_aln.bam'), 'w') as f:
        subprocess.call(['samtools', 'view', '-@', str(args.threads), '-bS', str(trinity_out + '/transcripts_aln.sam')], stdout=f, stderr=subprocess.DEVNULL)
    f.close()
    os.remove(str(trinity_out + '/transcripts_aln.sam'))

    with open(str(trinity_out + '/transcripts_aln.bed'), 'w') as f:
        subprocess.call(['bedtools', 'bamtobed', '-i', str(trinity_out + '/transcripts_aln.bam')], stdout=f)
    f.close()
    print(colored("Done!", "green", attrs=['bold']))


    clock = time()
    print(str(clock) + '\t' + "Masking transcripts with \"" + str(args.ref_TEs) + "\" database...")
    format = '.fa'
    test_fasta = re.search(format, str(args.ref_TEs))
    if test_fasta is None:
        subprocess.call(['RepeatMasker', str(trinity_out + '/Trinity.fasta'), '-species', str(args.ref_TEs), '-cutoff', '225', '-nolow', '-norna', '-s', '-par', str(args.threads)], stdout = subprocess.DEVNULL)
    else:
        subprocess.call(['RepeatMasker', str(trinity_out + '/Trinity.fasta'), '-lib', str(args.ref_TEs), '-cutoff', '225', '-nolow', '-norna', '-s', '-par', str(args.threads)], stdout = subprocess.DEVNULL)
    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(str(clock) + '\t' + "Identifying chimeric reads...")

    with open(str(trinity_out + '/Trinity.fasta.out'), 'r') as f:
        trinity_masked = f.read()
        trinity_masked = re.sub('^ +', '', trinity_masked, flags=re.M)
        trinity_masked = re.sub(' +', '\t', trinity_masked.strip())
    f.close()

    trinity_masked = StringIO(trinity_masked)
    trinity_masked = pd.read_table(trinity_masked, header=None, sep="\t", usecols=[4,5,6,9,2,8], names=['dot', 'transcript', 'start', 'end', 'strand', 'TE_fam']).iloc[2:]
    trinity_masked = trinity_masked[['transcript', 'start', 'end', 'TE_fam', 'dot', 'strand']]
    trinity_masked['strand'] = trinity_masked['strand'].replace('C', '-')
    trinity_masked['start'] = trinity_masked['start'].astype(int)
    trinity_masked['end'] = trinity_masked['end'].astype(int)
    trinity_masked['length'] = trinity_masked.apply(lambda x: x['end'] - x['start'], axis=1)
    trinity_masked = trinity_masked[trinity_masked.length > int(args.TE_length)]
    trinity_masked = trinity_masked[['transcript', 'start', 'end', 'TE_fam', 'dot', 'strand']]
    trinity_masked.to_csv(str(trinity_out + '/trinity_TEs.bed'), header = False, index = False, sep = '\t')

    intersection(str(trinity_out + '/transcripts_aln.bed'), str(trinity_out + '/trinity_TEs.bed'), str(trinity_out + '/overlapped.mbed'), str(args.overlap))
    print(colored("Done!", "green", attrs=['bold']))

    TE_reads = pd.read_csv(str(trinity_out + '/overlapped.mbed'), header=None, sep='\t', usecols=[3], names=['read_IDs']).replace('/1', '', regex=True).replace('/2', '', regex=True).value_counts().rename_axis('reads').reset_index(name='counts')
    TE_reads = TE_reads[TE_reads["counts"] == 1]#.to_csv(str(trinity_out + 'all_candidates'))
    TE_reads_list = TE_reads["reads"].to_list()

    overlapped = pd.read_csv(str(trinity_out + '/overlapped.mbed'), header=None, sep='\t',usecols=[0,3,9], names=['transcript_id','read_IDs','TE_fam'])
    overlapped[['read','mate']] = overlapped.read_IDs.str.split("/",expand=True)
    overlapped['mate'] = overlapped['mate'].replace('^', '/', regex=True)
    overlapped = overlapped[overlapped['read'].isin(TE_reads_list)]
    overlapped[['transcript_id', 'read_IDs', 'TE_fam']].to_csv(str(trinity_out + '/all_candidates.tsv'), header=None, index=False, sep="\t")

def singleton_crossing():
    from __main__ import aln_dir
    from __main__ import trinity_out
    from __main__ import group
    from __main__ import tmp

    isoforms = pd.read_csv(str(trinity_out + '/all_candidates.tsv'), header=None, sep='\t', usecols=[0], names=['transcript_id']).drop_duplicates()
    overlapped = pd.read_csv(str(trinity_out + '/all_candidates.tsv'), sep='\t', header=None, usecols=[0,1,2], names=['transcript_id','read_IDs','TE_fam'])

    transcript_reads = pd.read_csv(str(trinity_out + '/transcripts_aln.bed'), header=None, sep='\t',usecols=[0,3], names=['transcript_id','read_IDs'])
    transcript_reads[['read','mate']] = transcript_reads.read_IDs.str.split("/",expand=True)
    transcript_reads['mate'] = transcript_reads['mate'].replace('^', '/', regex=True)
    transcript_reads = transcript_reads[['transcript_id', 'read_IDs']]

    for row in isoforms.itertuples():
        chimeras = ' '

        transc_id = row[1]
        reads = overlapped.query('transcript_id == @transc_id').loc[:, ['read_IDs']].drop_duplicates()

        for read_id in reads.itertuples():
            ID = read_id[1]
            if '/1' in str(ID):
                mate = re.sub('/1','/2', ID)
            else:
                mate = re.sub('/2','/1', ID)
                chim_isoform = transcript_reads.query('read_IDs == @mate').query('transcript_id == @transc_id').loc[:, ['transcript_id']].drop_duplicates()
                if not chim_isoform.empty and chim_isoform is not None:
                    chim_isoform = chim_isoform.to_csv(header=None, index=False)

                    TE_fam = overlapped.query('read_IDs == @ID').loc[:, ['TE_fam']].drop_duplicates().to_csv(header=None, index=False, sep="\t")
                    if TE_fam is not None:
                        TE_fam = re.sub('\n', '', TE_fam); chim_isoform = re.sub('\n', '', chim_isoform)
                        if chim_isoform is not None:
                            isoform_TE = (str(chim_isoform) + '\t' + str(TE_fam))
                            chimeras += isoform_TE + '\n'

        if chimeras != ' ':
            chimeras = re.sub(' ', '',chimeras)
            chimeras_sIO = StringIO(chimeras)
            chimeras_df = pd.read_table(chimeras_sIO, header=None, sep='\t',usecols=[0,1], names=['transcript_id','TE_fam'])
            isoform_trinity = chimeras_df['transcript_id'].iloc[0]

            TE_freq = chimeras_df['TE_fam'].value_counts().rename_axis('reads').reset_index(name='counts').sort_values(by=['counts'], ascending=False).head(1).to_csv(header=None, index=False, sep="\t")
            TE_freq = re.sub('\n', '',TE_freq)
            with open(str(trinity_out + '/putative_chimeras.tsv'), 'a') as f:
                print(isoform_trinity + '\t' + TE_freq, file = f)
            f.close()

    clock = time()
    print(str(clock) + '\t' + "Performing blast to identify transcripts...")

    put_chim_isoforms = pd.read_table(str(trinity_out + '/putative_chimeras.tsv'), header=None, sep='\t',usecols=[0], names=['isoforms']).to_csv(str(trinity_out + '/putative_isoforms.lst'), header=None, index=False, sep="\t")
    with open(str(trinity_out + '/chim_trinity.fa'), 'w') as f:
        subprocess.call(['seqtk', 'subseq', str(trinity_out + '/Trinity.fasta'), str(trinity_out + '/putative_isoforms.lst')], stdout = f)
    f.close()

    subprocess.call(['makeblastdb', '-in', str(args.transcripts), '-dbtype', 'nucl', '-out', str(trinity_out + '/transcripts_db')], stdout=subprocess.DEVNULL)
    with open(str(trinity_out + '/blast-result.tsv'), 'w') as f:
        subprocess.call(['blastn', '-query', str(trinity_out + '/chim_trinity.fa'), '-db', str(trinity_out + '/transcripts_db'), '-outfmt', '6 qseqid sseqid length pident gaps mismatch qlen slen qstart qend sstart send evalue bitscore', '-num_threads', str(args.threads)], stdout =f)
    f.close()
    print(str(clock) + '\t' + colored("Done!", "green", attrs=['bold']))

    if check_file(str(trinity_out + '/blast-result.tsv')) == True:
        clock = time()
        print(str(clock) + '\t' + "Recovering the best matches...")
        blast_result_df = pd.read_table(str(trinity_out + '/blast-result.tsv'), header=None, sep='\t', names=['qseqid','sseqid','length','pident','gaps','mismatch','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'])

        isoform_IDs = blast_result_df.iloc[:, 0].drop_duplicates().to_list()
        for line in isoform_IDs:
            isoform_length = blast_result_df.query('qseqid == @line').loc[:, ['qlen']].head(1).to_csv(header=None, index=False, sep='\t')

            best_hit = blast_result_df.query('qseqid == @line').sort_values(by=['bitscore'], ascending=False).head(1)
            best_hit_ID = best_hit['sseqid'].to_csv(header=None, index=False, sep='\n').replace('\n', '')
            best_hit['length'] = best_hit['send'] - best_hit['sstart']
            best_hit_length = best_hit['length'].to_csv(header=None, index=False, sep='\n').replace('\n', '')

            identity = best_hit['pident'].to_csv(header=None, index=False, sep='\n').replace('\n', '')
            ref_len = best_hit['slen'].to_csv(header=None, index=False, sep='\n').replace('\n', '')
            chim_len = best_hit['qlen'].to_csv(header=None, index=False, sep='\n').replace('\n', '')

            transcript_gene = re.sub('_', '\t', best_hit_ID)
            perc = (int(best_hit_length) * 100) / int(isoform_length)
            if perc > int(args.identity):
                with open(str(trinity_out + '/IDs_isoforms.lst'), 'a') as isoforms:
                    print(str(line), file=isoforms)
                isoforms.close()
                with open(str(trinity_out + '/IDs_genes.lst'), 'a') as transcripts:
                    print(str(best_hit_ID), file = transcripts)
                transcripts.close()
                with open(str(trinity_out + '/blast_matches.tsv'), 'a') as output:
                    print(line + '\t' + transcript_gene + '\t' + identity + '\t' + chim_len + '\t' + ref_len + '\t' + best_hit_length, file=output)
                output.close()
        print(str(clock) + '\t' + colored("Done!", "green", attrs=['bold']))

        isoforms = pd.read_csv(str(trinity_out + '/blast_matches.tsv'), header=None, sep='\t', usecols=[0]).iloc[:, 0].drop_duplicates().to_list()
        put_chim_isoforms = pd.read_table(str(trinity_out + '/putative_chimeras.tsv'), header=None, sep='\t',usecols=[0,1,2], names=['isoform_id','TE_fam', 'cov'])
        blast_matches = pd.read_table(str(trinity_out + '/blast_matches.tsv'), header=None, sep='\t', names=['trinity_id','transcript_id','gene_id','pident','trinity_length','transcript_length','match_length'])

        for isoform in isoforms:
            TE_info = put_chim_isoforms.query('isoform_id == @isoform').loc[:, ['TE_fam', 'cov']].to_csv(header=None, index=False, sep='\t').replace('\n', '')
            blast_info = blast_matches.query('trinity_id == @isoform').to_csv(header=None, index=False, sep='\t').replace('\n', '')
            with open(str(trinity_out + '/' + group + '_transcriptome_evidence.tsv'), 'a') as output:
                print(blast_info + '\t' + TE_info, file=output)
            output.close()

        copy(str(trinity_out + '/' + group + '_transcriptome_evidence.tsv'), tmp)
    else:
        print(colored("No chimeric transcripts found by blast!", "red", attrs=['bold']))
        print(colored("You can try to decrease --identity parameter"))






#
