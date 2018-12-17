def mergeAnnotateLabel(mutect_infile, varscan_infile, tableannovar, varfiltpath, pair, humandb, buildver, buildver_mt, protocol_genome, protocol_mt, operations_genome, operations_mt, rsID_dict, processpath, labelled):

    import os
    import pandas as pd
    import subprocess as sp
    import json
    import re


    # create log folder if it doesn't already exists
    varfiltpath_logs = varfiltpath + 'logs/'
    os.makedirs(varfiltpath_logs, exist_ok=True)

    # load tables
    mutect = pd.read_table(mutect_infile, sep='\t', header=0, dtype='str')
    varscan = pd.read_table(varscan_infile, sep='\t', header=0, dtype='str')

    # filter MUTECT variants
    mutect = mutect.loc[(mutect.covered == 'COVERED') & (mutect.judgement == 'KEEP'),]

    # add 'detection_method' column
    mutect['detection_method'] = ['mutect'] * len(mutect)
    varscan['detection_method'] = ['varscan'] * len(varscan)

    # set column 'variant' as index
    mutect.index = mutect.variant
    varscan.index = varscan.variant

    # identify shared and unique variants between mutect and varscan lists
    unique_mutect = [v for v in mutect.variant if v not in varscan.variant]
    common = [v for v in varscan.variant if v in mutect.variant]

    # if a variant in varscan is in 'common', change its 'detection_method' field to 'both'
    varscan.loc[common, 'detection_method'] = 'both'

    # remove common variants from mutect
    mutect = mutect.loc[unique_mutect, :]


    ### MERGE MUTECT WITH VARSCAN, keeping the following columns order
    columns = ['chrom', 'start', 'end', 'ref', 'alt', 'n_vaf', 'n_depth', 't_vaf',
                't_depth', 'somatic_status', 'somaticFilter', 'tumor_strand_bias',
                'detection_method', 'variant']
    merged = pd.concat([mutect, varscan], sort=False).filter(items=columns)


    #### RUN ANNOVAR ####

    # set variables
    tableannovar = tableannovar
    genome_anno_infile = varfiltpath + pair + '_genome_annovar.tsv'
    mt_anno_infile = varfiltpath + pair + '_mt_annovar.tsv'
    humandb = humandb
    buildver = buildver
    buildver_mt = buildver_mt
    outlabel_genome = varfiltpath_logs + pair + '_annovar_genome'
    outlabel_mt = varfiltpath_logs + pair + '_annovar_mt'
    protocol_genome = protocol_genome
    protocol_mt = protocol_mt
    operations_genome = operations_genome
    operations_mt = operations_mt
    log_genome = varfiltpath_logs + pair + '_annovar_genome_err.log'
    log_mt = varfiltpath_logs + pair + '_annovar_mt_err.log'
    genome_anno_outfile = outlabel_genome + '.hg19_multianno.txt'
    mt_anno_outfile = outlabel_mt + '.GRCh37_MT_multianno.txt'

    # write first 5 columns of GENOMIC and MT variants to tables, as input for ANNOVAR
    merged.loc[merged.chrom != 'MT',].iloc[:,:5].to_csv(genome_anno_infile, sep='\t', index=False, header=False)
    merged.loc[merged.chrom == 'MT',].iloc[:,:5].to_csv(mt_anno_infile, sep='\t', index=False, header=False)

    # ANNOVAR on genomic variants
    cmd = '%s %s %s -buildver %s -out %s -remove -protocol %s -operation %s 2> %s' % (tableannovar, genome_anno_infile, humandb, buildver, outlabel_genome, protocol_genome, operations_genome, log_genome)
    sp.run(cmd, shell=True)

    # ANNOVAR on MT variants
    cmd = '%s %s %s -buildver %s -out %s -remove -protocol %s -operation %s 2> %s' % (tableannovar, mt_anno_infile, humandb, buildver_mt, outlabel_mt, protocol_mt, operations_mt, log_mt)
    sp.run(cmd, shell=True)

    # load ANNOVAR output tables
    genome = pd.read_table(genome_anno_outfile, sep='\t', header=0, dtype='str')
    mt = pd.read_table(mt_anno_outfile, sep='\t', header=0, dtype='str')

    # create 'variant' field and set to index (needed to later merge with other columns)
    genome.index = ['_'.join(row[:5].tolist() + [pair]) for i, row in genome.iterrows()]
    mt.index = ['_'.join(row[:5].tolist() + [pair]) for i, row in mt.iterrows()]



    ##################################################
    #### KEEP ONLY SOME COLUMNS OF ANNOVAR OUTPUT ####

    # set columns to keep
    genome_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                   'ExonicFunc.refGene', 'AAChange.refGene', 'cytoBand', 'genomicSuperDups',
                   'esp6500siv2_all', 'snp138', '1000g2014oct_all', 'ExAC_nontcga_ALL',
                   'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred',
                   'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred','MutationAssessor_score',
                   'MutationAssessor_pred', 'CLNALLELEID', 'CLNDN', 'CLNSIG', 'cosmic70']

    mt_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene',
               'ExonicFunc.ensGene', 'AAChange.ensGene']

    # filter out the other columns
    genome = genome.filter(items=genome_cols)
    mt = mt.filter(items=mt_cols)



    #############################
    #### CHANGE COLUMN NAMES ####

    # set new column names
    genome_new_cols = ['chrom', 'start', 'end', 'ref', 'alt', 'exon', 'gene', 'var_type',
                       'aa_change', 'cytoband', 'segdups', 'esp', 'dbsnp_id', '1000g',
                       'exac', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score',
                       'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred',
                       'MutationAssessor_score', 'MutationAssessor_pred', 'CLNALLELEID', 'CLNDN',
                       'CLNSIG', 'cosmic70']

    mt_new_cols = ['chrom', 'start', 'end', 'ref', 'alt', 'exon', 'gene', 'var_type',
                   'aa_change']

    # rename
    genome.columns = genome_new_cols
    mt.columns = mt_new_cols

    # merge rows GENOMIC and MT annotated variants, keeping the desired column order
    merged_annovar = pd.concat([genome, mt], sort=False).filter(items=genome_new_cols)

    # get final list of columns names
    columns = genome_new_cols + ['n_vaf', 'n_depth', 't_vaf', 't_depth', 'somatic_status',
                                 'somaticFilter', 'tumor_strand_bias', 'detection_method',
                                 'variant']

    # merge columns of every annotated variant with previous table (excluding first 5 columns)
    variants = pd.concat([merged_annovar, merged.iloc[:,5:]], axis=1, sort=False).filter(items=columns)



    #############################
    #### LABEL POLYMORPHISMS ####

    # load dictionary of polymorphic rsIDs
    with open(rsID_dict) as infile:
        rsID_dict = json.load(infile)


    # change to 'float' data in 3 columns
    variants = variants.astype({'1000g': 'float', 'esp': 'float', 'exac': 'float'})

    # set empty lists
    polymorphism = []

    # loop over rows
    for i, row in variants.iterrows():

        # set variables as 'no'
        polym = 'no'


        ## determine whether variants are polymorphisms or not ##

        # check if MAF in 1000genome, Exac, or Esp is >= 0.01
        if (row['1000g'] >= 0.01) | (row.esp >= 0.01) | (row.exac >= 0.01):
            polym = 'yes'

        # slice 8-mer prefix of dbsnp rsID
        prefix = str(row.dbsnp_id)[:8]

        # if the prefix is among the keys of the polymorphic rsIDs dictionary
        if prefix in rsID_dict.keys():
            # check if the current rsID in among the list of polymorphic rsIDs
            if row.dbsnp_id in rsID_dict[prefix]:
                polym = 'yes'

        # append variables to lists
        polymorphism.append(polym)

    # assign lists as new columns
    variants['polymorphism'] = polymorphism

    # write LABELLED variants to file
    variants.to_csv(labelled, sep='\t', header=True, index=False)


    # remove intermediate files
    sp.run('rm %s%s_*annovar*' % (varfiltpath, pair), shell=True)
    sp.run('rm %slogs/%s_*multianno*' % (varfiltpath, pair), shell=True)
    sp.run('rm %slogs/%s_*refGene*' % (varfiltpath, pair), shell=True)



mergeAnnotateLabel(snakemake.input['mutect_infile'],
                         snakemake.input['varscan_infile'],
                         snakemake.params['tableannovar'],
                         snakemake.params['varfiltpath'],
                         snakemake.params['pair'],
                         snakemake.params['humandb'],
                         snakemake.params['buildver'],
                         snakemake.params['buildver_mt'],
                         snakemake.params['protocol_genome'],
                         snakemake.params['protocol_mt'],
                         snakemake.params['operations_genome'],
                         snakemake.params['operations_mt'],
                         snakemake.params['rsID_dict'],
                         snakemake.params['processpath'],
                         snakemake.output['labelled'])
