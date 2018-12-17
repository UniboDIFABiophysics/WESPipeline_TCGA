def runVarscan(varscan, snv, indel, normal, tumor, min_avg_qual, strand_filter, min_var_freq, somatic_p_value, dataset_file, log):

    import subprocess as sp
    import os
    import pandas as pd

    # create fake .mpileup in case the real mpileup is empty
    # (only for MT variant call)
    if 'MT' in normal:

        fake_mpileup = '\t'.join(['MT', '1', 'A', '1', '.', '1'])

        if os.stat(normal).st_size == 0:
            with open(normal, 'w') as outfile:
                outfile.write(fake_mpileup)

        if os.stat(tumor).st_size == 0:
            with open(tumor, 'w') as outfile:
                outfile.write(fake_mpileup)


    # extract sample ids from paths
    tumor_id = tumor.split('/')[-1].split('.')[0]
    normal_id = normal.split('/')[-1].split('.')[0]

    # load dataset
    dataset = pd.read_table(dataset_file, sep = '\t',header=0, dtype='str', index_col='sample_id')

    # get purity values from dataset
    tumor_purity = dataset.loc[tumor_id, 'purity']
    normal_purity = dataset.loc[normal_id, 'purity']

    # if values are unknown, set standard ones
    if pd.isnull(tumor_purity): tumor_purity = '1'
    if pd.isnull(normal_purity): normal_purity = '0.8'

    # run VARSCAN
    cmd = 'java -jar %s somatic %s %s --output-snp %s --output-indel %s --min-avg-qual %s --strand_filter %s --min-var-freq %s --somatic-p-value %s --tumor-purity %s --normal-purity %s 2> %s' % (varscan, normal, tumor, snv, indel, min_avg_qual, strand_filter, min_var_freq, somatic_p_value, tumor_purity, normal_purity, log)
    sp.run(cmd, shell=True)

runVarscan(snakemake.input['varscan'],
                    snakemake.output['snv'],
                    snakemake.output['indel'],
                    snakemake.params['normal'],
                    snakemake.params['tumor'],
                    snakemake.params['min_avg_qual'],
                    snakemake.params['strand_filter'],
                    snakemake.params['min_var_freq'],
                    snakemake.params['somatic_p_value'],
                    snakemake.params['dataset_file'],
                    snakemake.params['log'])
