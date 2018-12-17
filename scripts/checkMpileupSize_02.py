def checkMpileupSize_02(varscan, normal, tumor, snv, indel, min_avg_qual, strand_filter, min_var_freq, somatic_p_value, log):

    import subprocess as sp
    import os

    if os.stat(normal).st_size == 0:
        with open(normal, 'w') as outfile:
            outfile.write('MT\t1\tA\t1\t.\t1')

    if os.stat(tumor).st_size == 0:
        with open(tumor, 'w') as outfile:
            outfile.write('MT\t1\tA\t1\t.\t1')

    cmd = 'java -jar %s somatic %s %s --output-snp %s --output-indel %s --min-avg-qual %s --strand_filter %s --min-var-freq %s --somatic-p-value %s 2> %s' % (varscan, normal, tumor, snv, indel, min_avg_qual, strand_filter, min_var_freq, somatic_p_value, log)
    sp.run(cmd, shell=True)

checkMpileupSize_02(snakemake.input['varscan'],
                    snakemake.input['normal'],
                    snakemake.input['tumor'],
                    snakemake.output['snv'],
                    snakemake.output['indel'],
                    snakemake.params['min_avg_qual'],
                    snakemake.params['strand_filter'],
                    snakemake.params['min_var_freq'],
                    snakemake.params['somatic_p_value'],
                    snakemake.params['log'])
