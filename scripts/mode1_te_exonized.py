from multiprocessing import Pool
import concurrent.futures
import shutil
from __main__ import *

#### Embedded functions
def prep_embedded(aln_dir):
    print("\n##########################\n## TE-exonized analysis ##\n##########################\n")

    ##Identify TEs within genes
    chim_exons = pd.read_csv(f'{aln_dir}/chim_exons.bed', sep = "\t", names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])

    ## Define list of gene IDs with overlapping with TEs at any region
    genes_TE_inside = intersection_any_bp(f'{aln_dir}/genes_total_expressed.bed', f'{aln_dir}/chimeric_TEs.bed')
    genes_TE_inside = genes_TE_inside["Name"].drop_duplicates().tolist()
    genes_TE_inside = set(genes_TE_inside)

    ## Select only chimeric exons from genes with TEs inside (remove exons from TE-init and TE-term)
    exons_with_TE_inside = chim_exons[chim_exons['Name'].isin(genes_TE_inside)]
    exons_with_TE_inside.to_csv(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed', sep = "\t", index = False, header = None)

    if os.stat(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed').st_size > 0:
        ## Create bed of exons with TEs embedded into
        embedded_exons_mbed = overlap_fraction_intersection_pyranges(f'{aln_dir}/chimeric_TEs.bed', f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed', fraction=1.0)
        if isinstance(embedded_exons_mbed, pd.DataFrame):
            embedded_exons_bed = embedded_exons_mbed[['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']]

            ## Create a list of gene IDs containing exons with emb. TEs
            genes_TE_embedded = embedded_exons_bed["Name_b"].drop_duplicates().tolist() #List of genes with embedded TEs on exons
            embedded_exons_bed.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']

            return genes_TE_embedded, embedded_exons_mbed, embedded_exons_bed
        else:
            return None

def intersector_fwd_rev_bed(gene_id, reads_gene_bed, df_emb_exons, TE_bed, reads_TEs_bed, overlap):
    chim_reads_exons = overlap_fraction_intersection_pyranges(reads_gene_bed, df_emb_exons, overlap)

    if isinstance(chim_reads_exons, pd.DataFrame):
        chim_reads_exons = get_IDs_from_bed(chim_reads_exons)
        chim_reads_exons = chim_reads_exons['Name'].str.replace(r'/[12]$', '', regex=True).drop_duplicates()
        
        chim_reads_TE = overlap_fraction_intersection_pyranges(reads_TEs_bed, TE_bed, overlap)
        if isinstance(chim_reads_TE, pd.DataFrame):
            ### Extract chimeric reads from TE
            chim_reads_TE = get_IDs_from_bed(chim_reads_TE)
            # Collapse /1 and /2 in-place
            chim_reads_TE['Name'] = chim_reads_TE['Name'].str.replace(r'/[12]$', '', regex=True)
            # Find duplicated names (True if appears >1 time)
            duplicated_mask = chim_reads_TE['Name'].duplicated(keep=False)
            # Keep only singleton rows (not duplicated)
            chim_reads_TE = chim_reads_TE.loc[~duplicated_mask, ['Name']]

            if not chim_reads_TE.empty:
                return chim_reads_exons, chim_reads_TE
            else:
                return None, None
        else:
            return None, None
    else:
        return None, None

def embedded_mp(gene_id, embedded_exons_bed, overlap_reads):
    from __main__ import aln_dir

    ### Extract all rows (TEs itersecting promoters)
    # exons_specific_gene = embedded_exons_bed.query('Name_b == @gene_id')[['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()    
    mask = embedded_exons_bed["Name_b"] == gene_id
    exons_specific_gene = embedded_exons_bed.loc[mask, ['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()
    
    ### Extract only strand
    strand_exon = exons_specific_gene['Strand_b'].values[0]

    ### Extract TE bed file
    maskTE = embedded_exons_bed["Name_b"] == gene_id
    TEs_specific_gene = embedded_exons_bed.loc[maskTE, ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']].drop_duplicates()
    
    
    if not TEs_specific_gene.empty:
        for row in TEs_specific_gene.itertuples():
            TE_coord_bed = pd.DataFrame([row._asdict()]).drop(columns='Index', errors='ignore')

            ### Check if embedded exons have chimeric reads
            if strand_exon == "+":
                reads_gene, reads_TE = intersector_fwd_rev_bed(gene_id, f'{aln_dir}/gene_reads_fwd.bed', exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_fwd.bed', overlap_reads)
            elif strand_exon == "-":
                reads_gene, reads_TE = intersector_fwd_rev_bed(gene_id, f'{aln_dir}/gene_reads_rev.bed', exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_rev.bed', overlap_reads)
                
            if isinstance(reads_gene, pd.DataFrame) and not reads_gene.empty \
                and isinstance(reads_TE, pd.DataFrame) and not reads_TE.empty:
                cov = reads_gene["Name"].isin(reads_TE["Name"]).sum()
                if cov > 0:
                    chr_TE = row[1]
                    s_TE = row[2]
                    e_TE = row[3]
                    TE_family = row[4]
                    TE_strd = row[6]
                
                    chr_gene   = exons_specific_gene['Chromosome'].values[0]
                    s_gene   = exons_specific_gene['Start_b'].values[0]
                    e_gene   = exons_specific_gene['End_b'].values[0]

                    with open(f'{aln_dir}/TE-exonized_embedded.tsv', 'a') as te_exo:
                        te_exo.write(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}\tEmbedded\n')



#### Overlapped functions      
def prep_overlapped(aln_dir, embedded_exons):
    ##Identify TEs overlapping exons
    if os.stat(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed').st_size > 0:
        exons_wsome_overlapping = pd.read_csv(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed', sep ='\t', usecols=[0,1,2,3,4,5],names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']).drop_duplicates() #All exons from genes with TEs >= 10% inside
    chimTEs = pd.read_csv(f'{aln_dir}/chimeric_TEs.bed', sep = "\t", names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])

    chimTEs_embedded = embedded_exons[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]
    chimTEs_embedded = chimTEs_embedded.drop_duplicates()
    chim_exons_embedded = embedded_exons[['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']]
    chim_exons_embedded.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']

    ### Subtract df with all exons with some overlap from the embedded exons df 
    overlapped_exons = exons_wsome_overlapping.merge(chim_exons_embedded.drop_duplicates(), how='left', indicator=True)
    overlapped_exons = overlapped_exons[overlapped_exons['_merge'] == 'left_only'].drop(columns=['_merge'])
    

    ### Subtract df with all chim TEs from the embedded TEs
    overlapped_TEs = chimTEs.merge(chimTEs_embedded.drop_duplicates(), how='left', indicator=True)
    overlapped_TEs = overlapped_TEs[overlapped_TEs['_merge'] == 'left_only'].drop(columns=['_merge'])

    overlap_TEs_exons_mbed = overlap_fraction_intersection_pyranges(overlapped_TEs, overlapped_exons, fraction=0.0)
    overlapped_gene_ids = overlap_TEs_exons_mbed["Name_b"].drop_duplicates().tolist()
    return overlap_TEs_exons_mbed, overlapped_gene_ids

def intersector_fwd_rev_bed_overlap(gene_id, reads_gene_bed, df_emb_exons, TE_coord_bed, reads_TEs_bed, overlap):
    chim_reads_exons = overlap_fraction_intersection_pyranges(reads_gene_bed, df_emb_exons, overlap)

    chim_reads_TE = None
    if isinstance(chim_reads_exons, pd.DataFrame) and not chim_reads_exons.empty:
        chim_reads_exons = get_IDs_from_bed(chim_reads_exons)
        chim_reads_exons = chim_reads_exons['Name'].str.replace(r'/[12]$', '', regex=True).drop_duplicates()
        
        if isinstance(TE_coord_bed, pd.DataFrame) and not TE_coord_bed.empty:
            chim_reads_TE = overlap_fraction_intersection_pyranges(reads_TEs_bed, TE_coord_bed, overlap)
            if isinstance(chim_reads_TE, pd.DataFrame):
                ### Extract chimeric reads from TE
                chim_reads_TE = get_IDs_from_bed(chim_reads_TE)
                # Collapse /1 and /2 in-place
                chim_reads_TE['Name'] = chim_reads_TE['Name'].str.replace(r'/[12]$', '', regex=True)
                # Find duplicated names (True if appears >1 time)
                duplicated_mask = chim_reads_TE['Name'].duplicated(keep=False)
                # Keep only singleton rows (not duplicated)
                chim_reads_TE = chim_reads_TE.loc[~duplicated_mask, ['Name']]
            else:
                return None, None
    if isinstance(chim_reads_TE, pd.DataFrame) and not chim_reads_TE.empty:
        return chim_reads_exons, chim_reads_TE
    else:
        return None, None

def overlapped_mp(gene_id, overlap_mbed, overlap_reads, genes_fwd, genes_rev):
    from __main__ import aln_dir
    
    ### Extract all rows (TEs itersecting promoters)
    mask = overlap_mbed["Name_b"] == gene_id
    exons_specific_gene = overlap_mbed.loc[mask, ['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()
    
    ### Extract only strand
    strand_exon = exons_specific_gene['Strand_b'].values[0]

    ### Extract TE bed file
    maskTE = overlap_mbed["Name_b"] == gene_id
    TEs_specific_gene = overlap_mbed.loc[maskTE, ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']].drop_duplicates()
    
    if not TEs_specific_gene.empty:
        for row in TEs_specific_gene.itertuples():
            TE_coord_bed = pd.DataFrame([row._asdict()]).drop(columns='Index', errors='ignore')

            # ### Check if embedded exons have chimeric reads
            if strand_exon == "+":
                reads_gene, reads_TE = intersector_fwd_rev_bed_overlap(gene_id, genes_fwd, exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_fwd.bed', overlap_reads)
            elif strand_exon == "-":
                reads_gene, reads_TE = intersector_fwd_rev_bed_overlap(gene_id, genes_rev, exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_rev.bed', overlap_reads)
            
            if isinstance(reads_gene, pd.DataFrame) and not reads_gene.empty \
                and isinstance(reads_TE, pd.DataFrame) and not reads_TE.empty:
                cov = reads_gene["Name"].isin(reads_TE["Name"]).sum()
                if cov > 0:
                    chr_TE = row[1]
                    s_TE = row[2]
                    e_TE = row[3]
                    TE_family = row[4]
                    TE_strd = row[6]
                
                    chr_gene   = exons_specific_gene['Chromosome'].values[0]
                    s_gene   = exons_specific_gene['Start_b'].values[0]
                    e_gene   = exons_specific_gene['End_b'].values[0]

                    with open(f'{aln_dir}/TE-exonized_overlapped.tsv', 'a') as te_exo:
                        te_exo.write(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}\tOverlapped\n')




#### Intronic functions
def prep_intronic(aln_dir):
    from __main__ import aln_dir

    ## Chimeric TEs vs expressed genes = Identify genes with TEs within
    genes_TE_inside_mbed = overlap_fraction_intersection_pyranges(f'{aln_dir}/chimeric_TEs.bed', f'{aln_dir}/genes_total_expressed.bed', fraction = 1.0)

    ## Subset TEs within genes
    TEs_within_genes = genes_TE_inside_mbed[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]

    ## Chimeric TEs vs all exons
    if os.stat(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed').st_size > 0:
        all_exons_genesTEinside = pd.read_csv(f'{aln_dir}/gene_TE_INSIDE-ALL_exons.bed', sep ='\t', usecols=[0,1,2,3,4,5],names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']).drop_duplicates() #All exons from genes with TEs >= 10% inside
    
    exons_overlap_TEs_mbed = overlap_fraction_intersection_pyranges(f'{aln_dir}/chimeric_TEs.bed', all_exons_genesTEinside, fraction = 0.0)
    TEs_some_overlap = exons_overlap_TEs_mbed[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]

    ### Subtract df with all TEs within genes from TEs with some overlap with exons 
    TEs_intron = TEs_within_genes.merge(TEs_some_overlap.drop_duplicates(), how='left', indicator=True)
    TEs_intron = TEs_intron[TEs_intron['_merge'] == 'left_only'].drop(columns=['_merge'])

    # TEs_intron.to_csv(f"{aln_dir}/TEs_intron.tsv", sep = "\t", index = False)

    ### All exons from genes with chimeric TEs at intron
    genes_TE_intron = overlap_fraction_intersection_pyranges(f'{aln_dir}/genes_total_expressed.bed', TEs_intron, fraction = 0.0)
    # print("genes_TE_intron")
    # print(genes_TE_intron)
    # print("genes_TE_inside_mbed")
    # print(genes_TE_inside_mbed)
    geneIDs_TE_intron = get_IDs_from_bed(genes_TE_intron)
    exons_genes_TEs_intron = all_exons_genesTEinside[all_exons_genesTEinside["Name"].isin(geneIDs_TE_intron["Name"])]
    list_geneIDs_TE_intron = geneIDs_TE_intron["Name"].drop_duplicates().to_list()

    return list_geneIDs_TE_intron, exons_genes_TEs_intron, genes_TE_intron

def intersector_fwd_rev_bed_intron(gene_id, reads_gene_bed, df_emb_exons, TE_coord_bed, reads_TEs_bed, overlap):
    chim_reads_exons = overlap_fraction_intersection_pyranges(reads_gene_bed, df_emb_exons, overlap)
    chim_reads_TE = None
    if isinstance(chim_reads_exons, pd.DataFrame):
        chim_reads_exons = get_IDs_from_bed(chim_reads_exons)
        chim_reads_exons = chim_reads_exons[['Name']].replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
       
        chim_reads_TE = overlap_fraction_intersection_pyranges(reads_TEs_bed, TE_coord_bed, overlap)
        if isinstance(chim_reads_TE, pd.DataFrame):
            ### Extract chimeric reads from TE
            chim_reads_TE = get_IDs_from_bed(chim_reads_TE)
            # Collapse /1 and /2 in-place
            chim_reads_TE['Name'] = chim_reads_TE['Name'].str.replace(r'/[12]$', '', regex=True)
            # Find duplicated names (True if appears >1 time)
            duplicated_mask = chim_reads_TE['Name'].duplicated(keep=False)
            # Keep only singleton rows (not duplicated)
            chim_reads_TE = chim_reads_TE.loc[~duplicated_mask, ['Name']]
        else:
            return None, None
    if isinstance(chim_reads_exons, pd.DataFrame):
        return chim_reads_exons, chim_reads_TE
    else:
        return None, None
    
def intronic_mp(gene_id, all_exons_bed, mbed_genes, overlap_reads):
    from __main__ import aln_dir

    ### Extract all rows (TEs itersecting promoters)
    mask = all_exons_bed["Name"] == gene_id
    exons_specific_gene = all_exons_bed.loc[mask, ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']].drop_duplicates()
    if isinstance(exons_specific_gene, pd.DataFrame) and not exons_specific_gene.empty:
        ## Extract only strand
        strand_exon = exons_specific_gene['Strand'].values[0]

        ### Extract TE bed file
        maskTE = mbed_genes["Name"] == gene_id
        TEs_specific_gene = mbed_genes.loc[maskTE, ['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()
    
        if not TEs_specific_gene.empty:
            for row in TEs_specific_gene.itertuples():
                TE_coord_bed = pd.DataFrame([row._asdict()]).drop(columns='Index', errors='ignore')

                ### Check if embedded exons have chimeric reads
                if strand_exon == "+":
                    reads_gene, reads_TE = intersector_fwd_rev_bed_intron(gene_id, f'{aln_dir}/gene_reads_fwd.bed', exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_fwd.bed', overlap_reads)
                elif strand_exon == "-":
                    reads_gene, reads_TE = intersector_fwd_rev_bed_intron(gene_id, f'{aln_dir}/gene_reads_rev.bed', exons_specific_gene, TE_coord_bed, f'{aln_dir}/TE_reads_rev.bed', overlap_reads)
                    
                if isinstance(reads_gene, pd.DataFrame) and not reads_gene.empty \
                and isinstance(reads_TE, pd.DataFrame) and not reads_TE.empty:
                    cov = reads_gene["Name"].isin(reads_TE["Name"]).sum()
                    if cov > 0:
                        chr_TE = row[1]
                        s_TE = row[2]
                        e_TE = row[3]
                        TE_family = row[4]
                        TE_strd = row[6]
                
                        chr_gene   = exons_specific_gene['Chromosome'].values[0]
                        s_gene   = exons_specific_gene['Start'].values[0]
                        e_gene   = exons_specific_gene['End'].values[0]

                        with open(f'{aln_dir}/TE-exonized_intronic.tsv', 'a') as te_exo:
                            te_exo.write(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}\n')



### Multi-processor funcions
def multicore_emb_exon(emb_ids, te_exon_bed, threads, overlap_reads):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import tmp

    ### TE-exonized EMBEDDED chimeras
    clock = time()
    print(f'{clock}\tSearching Embedded-TE chimeras...')

    if os.path.exists(f'{aln_dir}/TE-exonized_embedded.tsv'):
        os.remove(f'{aln_dir}/TE-exonized_embedded.tsv')
    
    ### Run multi-threads te-exonized embedded function
    if len(emb_ids) > 0:
        with Pool(threads) as pool:
            pool.starmap(embedded_mp, [(gene_id, te_exon_bed, overlap_reads) for gene_id in emb_ids])       
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(colored('There are no TEs embedded into exons with chimeric reads!', "red"))

    if os.path.exists(f'{aln_dir}/TE-exonized_embedded.tsv'):
        if os.stat(f'{aln_dir}/TE-exonized_embedded.tsv').st_size > 0:
            os.rename(f'{aln_dir}/TE-exonized_embedded.tsv', f'{aln_dir}/TE-exonized_embedded-{group}.tsv')
            shutil.move(f'{aln_dir}/TE-exonized_embedded-{group}.tsv', f'{tmp}/TE-exonized_embedded-{group}.tsv')

def multicore_overl_exon(overl_ids, overl_mbed, threads, overlap_reads):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import tmp

    ### TE-exonized OVERLAPPED chimeras
    clock = time()
    print(f'{clock}\tSearching Overlapped-TE chimeras...')

    if os.path.exists(f'{aln_dir}/TE-exonized_overlapped.tsv'):
        os.remove(f'{aln_dir}/TE-exonized_overlapped.tsv')

    ## Import gene reads 
    gene_fwd_reads = pd.read_csv(f'{aln_dir}/gene_reads_fwd.bed', sep = "\t", header=None)
    gene_rev_reads = pd.read_csv(f'{aln_dir}/gene_reads_rev.bed', sep = "\t", header=None)

    ### Run multi-threads te-exonized embedded function
    if len(overl_ids) > 0:
        with Pool(threads) as pool:
            pool.starmap(overlapped_mp, [(gene_id, overl_mbed, overlap_reads, gene_fwd_reads, gene_rev_reads) for gene_id in overl_ids])       
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(colored('There are no TEs embedded into exons with chimeric reads!', "red"))

    if os.path.exists(f'{aln_dir}/TE-exonized_overlapped.tsv'):
        if os.stat(f'{aln_dir}/TE-exonized_overlapped.tsv').st_size > 0:
            os.rename(f'{aln_dir}/TE-exonized_overlapped.tsv', f'{aln_dir}/TE-exonized_overlapped-{group}.tsv')
            shutil.move(f'{aln_dir}/TE-exonized_overlapped-{group}.tsv', f'{tmp}/TE-exonized_overlapped-{group}.tsv')







    # if os.path.exists(list) == True:
    #     with open(str(aln_dir + '/exon_overlapped.lst')) as f:
    #         all_overlap_IDs = f.read().splitlines()
    #     f.close
    #     pool = Pool(processes=int(args.threads))
    #     pool.map(overlapped_mp, all_overlap_IDs)
    #     pool.close()
    #     print(colored("Done!", "green", attrs=['bold']))
    #     pybedtools.cleanup(remove_all=True)
    # else:
    #     print(colored('There are no TEs overlapping exons with chimeric reads!', "red"))



    # ### TE-exonized INTRONIC chimeras
    # clock = time()
    # print(str(clock) + '\t' + "Searching Intronic-TE chimeras...")
    # list = str(aln_dir + '/exon_intron_list.lst')

    # if os.path.exists(list) == True:
    #     with open(str(aln_dir + '/exon_intron_list.lst')) as f:
    #         all_gene_IDs = f.read().splitlines()
    #     f.close
    #     pool = Pool(processes=int(args.threads))
    #     pool.map(intronic_mp, all_gene_IDs)
    #     pool.close()
    #     print(colored("Done!", "green", attrs=['bold']))
    #     pybedtools.cleanup(remove_all=True)
    # else:
    #     print(colored('There are no TEs into introns with chimeric reads!', "red"))

def multicore_intron_exon(intron_ids, exons_bed, mbed_genes, threads, overlap_reads):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import tmp

    ### TE-exonized OVERLAPPED chimeras
    clock = time()
    print(f'{clock}\tSearching Intronic-TE chimeras...')

    if os.path.exists(f'{aln_dir}/TE-exonized_intronic.tsv'):
        os.remove(f'{aln_dir}/TE-exonized_intronic.tsv')

    ### Run multi-threads te-exonized embedded function
    if len(intron_ids) > 0:
        with Pool(threads) as pool:
            pool.starmap(intronic_mp, [(gene_id, exons_bed, mbed_genes, overlap_reads) for gene_id in intron_ids])       
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(colored('There are no TEs embedded into exons with chimeric reads!', "red"))

    if os.path.exists(f'{aln_dir}/TE-exonized_intronic.tsv'):
        if os.stat(f'{aln_dir}/TE-exonized_intronic.tsv').st_size > 0:
            os.rename(f'{aln_dir}/TE-exonized_intronic.tsv', f'{aln_dir}/TE-exonized_intronic-{group}.tsv')
            shutil.move(f'{aln_dir}/TE-exonized_intronic-{group}.tsv', f'{tmp}/TE-exonized_intronic-{group}.tsv')