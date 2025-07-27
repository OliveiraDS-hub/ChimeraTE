from __main__ import *

def stranded_reads(strand_option, aln_dir, threads):
    if strand_option == "fwd-stranded":
        pair_sense = ["rev", "fwd"]
    elif strand_option == "rf-stranded":
        pair_sense = ["fwd", "rev"]

    for index, sense in enumerate(pair_sense):
        print(f"Creating bed file for {sense} reads...")
        if index == 0:
            ### First strand
            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '128', '-F', '16', f"{aln_dir}/accepted_hits.bam", '-o', f"{aln_dir}/{sense}1_f.bam"], stderr=subprocess.DEVNULL)

            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '80', f"{aln_dir}/accepted_hits.bam", '-o', f"{aln_dir}/{sense}2_f.bam"], stderr=subprocess.DEVNULL)
        elif index == 1:
            ### Second strand
            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '144', f"{aln_dir}/accepted_hits.bam", '-o', f"{aln_dir}/{sense}1_f.bam"], stderr=subprocess.DEVNULL)

            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '64', '-F', '16', f"{aln_dir}/accepted_hits.bam", '-o', f"{aln_dir}/{sense}2_f.bam"], stderr=subprocess.DEVNULL)

        ### Combine alignments that originate on the forward strand.
        subprocess.call(['samtools', 'merge', '-@', str(threads), '-f', f"{aln_dir}/{sense}.bam", f"{aln_dir}/{sense}1_f.bam", f"{aln_dir}/{sense}2_f.bam"])

        subprocess.run('bedtools bamtobed -split -i {}/{}.bam > {}/{}.bed'.format(aln_dir, sense, aln_dir, sense), shell = True)

def overlap_reads(df, read_ids, aln_dir=None, output=None):
        # Ensure fast lookup by using a set
        read_ids_set = set(read_ids)

        # Define col names
        df.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']

        # Normalize read names (remove /1 and /2)
        df['ID_mod'] = df["Name"].str.replace('/1', '', regex=False).str.replace('/2', '', regex=False)

        # Filter rows with matching IDs
        filtered = df[df['ID_mod'].isin(read_ids_set)]

        # Columns to keep
        cols = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']

        # Output
        if output is None:
            filtered.to_csv(sys.stdout, header=None, index=False, sep="\t", columns=cols)
        else:
            filtered.to_csv(f"{aln_dir}/{output}.bed", header=None, index=False, sep="\t", columns=cols)

def alignment_func(out_dir,group,aln_dir,mate1,mate2, threads, strandness, tmp_dir, fpkm, index=False):
    print(f'Running analysis with ------------------------------------------> {group}\n')
    clock = time()
    print(f'{clock}\tPerforming alignment')

    prop_cpu = int(int(threads) - (int(threads) * 0.2))
    
    ### STAR alignment | output = sorted bam file
    # Check write permission
    if not os.access(aln_dir, os.W_OK):
        raise PermissionError(f"No write permission for {aln_dir}")
    
    if os.path.exists(f'{out_dir}/{group}_Aligned.sortedByCoord.out.bam') and \
        os.path.getsize(f'{out_dir}/{group}_Aligned.sortedByCoord.out.bam') == 0:
        os.remove(f'{out_dir}/{group}_Aligned.sortedByCoord.out.bam')
    
    if not os.path.exists(f"{aln_dir}/{group}_Aligned.sortedByCoord.out.bam"):
        if not index == False:
            if mate1.endswith('gz'):
                subprocess.call(['STAR', '--genomeDir', str(index), '--runThreadN', str(prop_cpu), "--readFilesCommand", "zcat", \
                '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', "BAM", "SortedByCoordinate", '--outTmpDir', f"{aln_dir}/{group}_STARtmp", "--outFileNamePrefix", f"{aln_dir}/{group}_"], stdout=subprocess.DEVNULL)
            else:
                subprocess.call(['STAR', '--genomeDir', str(index), '--runThreadN', str(prop_cpu), \
                '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', "BAM", "SortedByCoordinate", '--outTmpDir', f"{aln_dir}/{group}_STARtmp", "--outFileNamePrefix", f"{aln_dir}/{group}_"], stdout=subprocess.DEVNULL)
        else:
            if mate1.endswith('gz'):
                subprocess.call(['STAR', '--genomeDir', f"{out_dir}/index", '--runThreadN', str(prop_cpu), "--readFilesCommand", "zcat", \
                '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', "BAM", "SortedByCoordinate", '--outTmpDir', f"{aln_dir}/{group}_STARtmp", "--outFileNamePrefix", f"{aln_dir}/{group}_"], stdout=subprocess.DEVNULL)
            else:
                subprocess.call(['STAR', '--genomeDir', f"{out_dir}/index", '--runThreadN', str(prop_cpu), \
                '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', "BAM", "SortedByCoordinate",'--outTmpDir', f"{aln_dir}/{group}_STARtmp", "--outFileNamePrefix", f"{aln_dir}/{group}_"], stdout=subprocess.DEVNULL)
    else:
        print(f'BAM file found: {aln_dir}/{group}_Aligned.sortedByCoord.out.bam\t Skipping alignment')

    if os.path.exists(f'{out_dir}/accepted_hits.bed'):
        os.remove(f'{out_dir}/accepted_hits.bed')

    if not os.path.exists(f'{aln_dir}/accepted_hits.bed'):
        ### samtools conversion | output = bam file with concordant paired reads
        with open(f"{aln_dir}/accepted_hits.bam", 'w') as bam_file:
            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-q', '255', f"{aln_dir}/{group}_Aligned.sortedByCoord.out.bam"], stdout=bam_file)
        bam_file.close
        ### samtools conversion | output = bed file with coordiantes of aligned reads
        with open(f'{aln_dir}/accepted_hits.bed', 'w') as bed_file:
            subprocess.call(['bedtools', 'bamtobed', '-split', '-i', f'{aln_dir}/accepted_hits.bam'], stdout=bed_file)
        bed_file.close

    ### Creating bed files per strand
    if not os.path.exists(f'{aln_dir}/rev.bed') and not os.path.exists(f'{aln_dir}/fwd.bed'):
        stranded_reads(strandness, aln_dir, threads)
        print(colored("Done!", "green", attrs=['bold']))
    




    clock = time()
    print(f'{clock}\tGenes expression')
    ### Calculate gene expression in fpkm with cufflinks | output = fpkm/genes.fpkm_tracking
    if not os.path.exists(f'{out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking'):
        fpkm_dir = create_dir(f'{out_dir}/{group}/alignment/fpkm')

        if str(strandness) == str("rf-stranded"):
            subprocess.call(['cufflinks', f'{aln_dir}/accepted_hits.bam', '-p', str(threads), '-G', f'{tmp_dir}/gtf_file.gtf', '-o', str(fpkm_dir), '--quiet', '--library-type', 'fr-firststrand'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        else:
            subprocess.call(['cufflinks', f'{aln_dir}/accepted_hits.bam', '-p', str(threads), '-G', f'{tmp_dir}/gtf_file.gtf', '-o', str(fpkm_dir), '--quiet', '--library-type', 'fr-secondstrand'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print(colored("Done!", "green", attrs=['bold']))

    ### Filter gene expression based on minimum fpkm | output = genes_expressed_IDs.lst
    if os.path.isfile(f'{out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking') and os.path.getsize(f'{out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking') > 0:
        fpkm_gene = pd.read_csv(f'{out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking', sep="\t", usecols=[0,11])
        fpkm_gene['tracking_id'] = fpkm_gene['tracking_id'].str.replace("gene-", '')
        list_genes = (fpkm_gene.loc[fpkm_gene['FPKM_conf_hi'] > int(fpkm)]
            .loc[:, ['tracking_id']]
            .drop_duplicates()['tracking_id']
            .tolist())
    else:
        raise RuntimeError(f"Expected output from cufflinks: {out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking was not found or is empty!"); return False
    
    if len(list_genes) > 0:
        try:
            clock = time()
            print(f'{clock}\tStrand-specific expression analysis')

            ### Create bed file only for expressed genes | output = genes_total_expressed.bed
            overlap(f'{tmp_dir}/gene_coord.bed', list_genes, f'{aln_dir}/genes_total_expressed.bed')

            ### Create bed file only for expressed TEs (at least 1 read) | output = expressed_TEs.bed
            if not os.path.exists(f'{aln_dir}/expressed_TEs.bed'):
                intersection_any_bp(f'{tmp_dir}/TE_file.bed', f'{aln_dir}/accepted_hits.bed', f'{aln_dir}/expressed_TEs.bed')

            ### Extract TE read IDs from expressed copies | output = hold reads into variables
            fwd_TE = intersection_any_bp(f'{aln_dir}/fwd.bed', f'{tmp_dir}/TE_file.bed', False)
            rev_TE = intersection_any_bp( f'{aln_dir}/rev.bed', f'{tmp_dir}/TE_file.bed', False)

            # ### Merge fwd and rev reads | output = list with read IDs
            merged_TE = pd.concat([fwd_TE, rev_TE])
            read_list_TE = merged_TE["Name"].str.replace('/1', '').str.replace('/2', '').drop_duplicates().to_list()
            
            # ### Extract Gene read IDs from expressed copies | output = hold reads into variables
            fwd_gene = intersection_any_bp(f'{aln_dir}/fwd.bed', f'{aln_dir}/genes_total_expressed.bed', False)
            rev_gene = intersection_any_bp(f'{aln_dir}/rev.bed', f'{aln_dir}/genes_total_expressed.bed', False)

            ### Merge fwd and rev reads | output = TE_reads.lst
            merged_gene = pd.concat([fwd_gene, rev_gene])
            # merged_gene["Name"].str.replace('/1', '').str.replace('/2', '').drop_duplicates().to_csv(aln_dir + '/gene_reads.lst', header=None, index=False)
            read_list_gene = merged_gene["Name"].str.replace('/1', '').str.replace('/2', '').drop_duplicates().to_list()


            print(f'total reads from TEs: {len(read_list_TE)}')
            print(f'total reads from genes: {len(read_list_gene)}')
            print(colored("Done!", "green", attrs=['bold']))
        except ValueError as e:
            print("Failed to intersect foward and reverse reads", e)


        if not os.path.exists(f'{aln_dir}/chim_exons.bed') and not os.path.exists(f'{aln_dir}/chimeric_TEs.bed'):
            clock = time()
            print(f'{clock}\tChimeric reads pairs identification')
            
            ### Load TE and gene reads as a sets
            te_reads_set = set(read_list_TE)
            gene_reads_set = set(read_list_gene)

            ### Find intersection
            chim_reads_set = te_reads_set & gene_reads_set

            ### Convert back to list if needed
            chim_reads = list(chim_reads_set)


            ### Create bed file with TEs' chimeric reads fwd and rev | output = fwd/rev.bed files
            overlap_reads(fwd_TE, chim_reads, aln_dir, "TE_reads_fwd")
            overlap_reads(rev_TE, chim_reads, aln_dir, "TE_reads_rev")

            ### Create bed file with genes' chimeric reads fwd and rev | output = fwd/rev.bed files
            overlap_reads(fwd_gene, chim_reads, aln_dir, "gene_reads_fwd")
            overlap_reads(rev_gene, chim_reads, aln_dir, "gene_reads_rev")


            ### Create bed file only for chimeric reads | output = chim_reads.bed
            accepted_hits = pd.read_csv(f'{aln_dir}/accepted_hits.bed', sep = "\t", names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
            overlap_reads(accepted_hits, chim_reads, aln_dir, 'chim_reads')

            ### Detect exons with chimeric reads | output = chim_exons.bed
            intersection_any_bp(f'{tmp_dir}/exon_file.bed', f'{aln_dir}/chim_reads.bed', f'{aln_dir}/chim_exons.bed') 
            
            ### Detect TEs with chimeric reads | output = chimeric_TEs.bed
            intersection_any_bp(f'{aln_dir}/expressed_TEs.bed', f'{aln_dir}/chim_reads.bed', f'{aln_dir}/chimeric_TEs.bed') 
            return True
            print(colored("Done!", "green", attrs=['bold']))
        else:
            print(f'chimeric reads for TEs and genes have been found!\tSkipping it...'); return True
    else:
        print(f"{group} has no genes with fpkm >= {fpkm}!\t Exiting...")
        return None
    