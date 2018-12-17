
#% TODO: dbsnp index???

shell.executable('/bin/bash')

import pandas as pd
import os
import subprocess as sp

###########################
##### set main paths ######

# get main path (home)
homepath = os.path.expanduser('~') + '/'

# CONDA
condapath = homepath + config['conda']

# SCRIPTS
scriptpath = config['scripts']


# build and create PROCESSPATH
processpath = homepath + config['process_path'] + config['outfolder'] + '/'


################################################################################

# load dataset file and set sample names as indeces
dataset_file = homepath + config['dataset']
dataset = pd.read_table(dataset_file, sep = '\t',header=0, dtype='str')
dataset = dataset.set_index('sample_id')

# References
rsID_dict = homepath + config['rsID_dict']

hg = homepath + config['fasta']['genome'] # Human Genome Reference
MT = homepath + config['fasta']['MT']
hg_indexes = [hg+'.bwt', hg+'.pac', hg+'.amb', hg+'.ann', hg+'.sa'] # Indexes for hg

dbsnp = homepath+ config['ref-files']['dbsnp'] # SNP database
cosmic = homepath+ config['ref-files']['cosmic'] # Catalog of somatic mutation in cancer
humandb = homepath+ config['ref-files']['humandb'] # Annovar human databases folder
buildver = config['ref-files']['buildver'] # Set build version
buildver_mt = config['ref-files']['buildver_mt'] # Set parameter command for annotating mitochondria variants


# Varscan arguments
min_avg_qual = config['varscan_arguments']['min_avg_qual']
strand_filter = config['varscan_arguments']['strand_filter']
min_var_freq = config['varscan_arguments']['min_var_freq']
somatic_p_value = config['varscan_arguments']['somatic_p_value']

# Varscan somaticFilter arguments
min_coverage = config['somaticFilter']['min_coverage']
min_reads2 = config['somaticFilter']['min_reads2']
min_var_freq_sf = config['somaticFilter']['min_var_freq']


# Softwares
gatk = homepath+ config['softwares']['gatk']
muTect = homepath+ config['softwares']['muTect']
annovar = homepath + config['folders']['annovar']
annotate = homepath+ config['softwares']['annotate']
tableannovar = homepath+ config['softwares']['tableannovar']
varscan = homepath + config['softwares']['varscan']

# bed files
nextexp_bed = homepath + config['bed']['nextera_expanded'] # Set target intervals for exome analysis
nextera_bed = homepath + config['bed']['nextera']
MT_bed = homepath + config['bed']['MT']
truseq_bed = homepath + config['bed']['truseq']
truseq_rapid_bed = homepath + config['bed']['truseq_rapid']
nextera_dna_exome_bed = homepath + config['bed']['nextera_dna_exome']


# Annovar databases
annovar_dbs = [homepath + config['annovar_dbs']['hg19_refGene'],
               homepath + config['annovar_dbs']['hg19_cytoBand'],
               homepath + config['annovar_dbs']['hg19_gSd'],
               homepath + config['annovar_dbs']['hg19_esp'],
               homepath + config['annovar_dbs']['hg19_snp138'],
               homepath + config['annovar_dbs']['hg19_1000g2014oct'],
               homepath + config['annovar_dbs']['hg19_exac03nontcga'],
               homepath + config['annovar_dbs']['hg19_ljb26_all'],
               homepath + config['annovar_dbs']['hg19_clinvar_20180603'],
               homepath + config['annovar_dbs']['hg19_cosmic70'],
               ]

protocol_genome = config['protocols']['genome']
protocol_mt = config['protocols']['mt']
operations_genome = config['operations']['genome']
operations_mt = config['operations']['mt']


# Label's parameters
n_sim = config['n_sim']
cpu_type = config['cpu_type']
thrs = config['threads']
n_cpu = config['n_cpu']


# move sample list file to processpath
input_list = homepath + config['input_list']
sp.run('cp %s %s' % (input_list, processpath), shell=True)

# load table
df_pair = pd.read_table(input_list, sep= '\t', header=0, dtype='str')


alarms = []
pairs = []
kits = []
# KIT CHECKPOINT
for i, row in df_pair.iterrows():

    pair = '_'.join([row.tumor, row.normal])

    # get WES library kit for tumor and normal
    kit_tumor = dataset.loc[row.tumor, 'kit_wes']
    kit_normal = dataset.loc[row.normal, 'kit_wes']

    if kit_tumor != kit_normal:
        alarms.append((pair, 'different_kit_wes'))

#    kit = '_'.join([kit_tumor, kit_normal])

    pairs.append(pair)
    kits.append(kit_tumor)

# complete pair DF
df_pair.index = pairs
df_pair['kit'] = kits

# create alarms DF and save to file
alarms = pd.DataFrame(alarms, columns=['pair', 'alarm'])
alarms.to_csv(processpath + 'alarms.tsv', sep='\t', index=False)

# create sample DF
df_sample = pd.DataFrame({'sample': df_pair.tumor.tolist() + df_pair.normal.tolist(),
                          'directory': df_pair.tumor_directory.tolist() + df_pair.normal_directory.tolist()})
df_sample = df_sample.set_index('sample')


# create list of every sample
samples = [s for p in pairs for s in p.split('_')]
samples = list(set(samples))
# Wildcard costrains necessary for search only certain names
wildcard_constraints:
    sample = "("+"|".join(samples)+")",
    pair = "("+"|".join(pairs)+")",



### set folder tree

alignpath_genome = processpath + '03_alignment_genome/'
alignpath_genomebqsr = alignpath_genome + '02_bqsr/'

alignpath_MT = processpath + '05_alignment_MT/'
alignpath_MTbqsr = alignpath_MT + '02_bqsr/'

mutectpath_genome = processpath + '06_mutect_genome/'
mutectpath_genome_logs = mutectpath_genome + 'logs/'
mutectpath_MT = processpath + '07_mutect_MT/'
mutectpath_MT_logs = mutectpath_MT + 'logs/'

varscanpath_genome = processpath + '08_varscan_genome/'
varscanpath_MT = processpath + '09_varscan_MT/'
mpileup_genome = varscanpath_genome + '01_mpileup/'
mpileup_MT = varscanpath_MT + '01_mpileup/'
varscan_genome_logs = varscanpath_genome + 'logs/'
varscan_MT_logs = varscanpath_MT + 'logs/'

variants_filter = processpath + '10_variants_filter/'

benchmarkpath = processpath + 'benchmarks/'





def get_bed_pair(wildcards):
    kit = df_pair.loc[wildcards, 'kit']
    bed = homepath + config['bed'][kit] + '_fixed.bed'
    return(bed)

def get_bam_bai_pair(wildcards, pairing, genome_type, extension):
    sample = df_pair.loc[wildcards, pairing]

    d = {'nuclear': '03_alignment_genome/02_bqsr/',
         'MT': '05_alignment_MT/02_bqsr/'}
    store_dir = pairing + '_directory'
    store_dir = df_pair.loc[wildcards, store_dir] + d[genome_type]
    store_dir = '%s%s%s' %(store_dir, sample, extension)

    d = {'nuclear': alignpath_genomebqsr,
         'MT': alignpath_MTbqsr}
    wes_dir = '%s%s/%s%s' % (d[genome_type], wildcards, sample, extension)

    d = {'store_dir': store_dir, 'wes_dir': wes_dir}
    return(d)

def get_bed(wildcards):
    kit = dataset.loc[wildcards, 'kit_wes']
    bed = homepath + config['bed'][kit] + '_fixed.bed'
    return(bed)

def get_mpileup(wildcards, pairing, genome_type):
    sample = df_pair.loc[wildcards, pairing]
    d = {'nuclear': mpileup_genome,
         'MT': mpileup_MT}
    mpileup = '%s%s/%s.mpileup' % (d[genome_type], wildcards, sample)
    return(mpileup)






##########################################
#           PIPELINE BEGINNING           #
##########################################


rule all:
    '''
      PIPELINE ENDING
    '''
    input:
        expand(mutectpath_genome + '{pair}_cov_wig.log.gz', pair=pairs),
        expand(mutectpath_MT + '{pair}_cov_wig.log.gz', pair=pairs),
        expand(variants_filter + '{pair}_labelled_variants.tsv', pair=pairs),
    message: ' THE END '
    run:
        pass


###########################################################################################



rule retrieveBAM_nuclear:
    '''
    Download the necessary nuclear BAM and BAI from the storepath
    '''
    input:
        normal_bam = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'nuclear', '.bam')['store_dir'],
        normal_bai = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'nuclear', '.bai')['store_dir'],
        tumor_bam = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'nuclear', '.bam')['store_dir'],
        tumor_bai = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'nuclear', '.bai')['store_dir'],
    output:
        temp(alignpath_genomebqsr + '{pair}/')
    resources:
        disk = 1
    benchmark:
        benchmarkpath + 'benchmark_retrieve_BAM_nuclear_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Dowloading nuclear BAM/BAI from storepath'
    shell:
        'rsync -a {input.normal_bam} {input.normal_bai} {input.tumor_bam} {input.tumor_bai} {output}'


rule retrieveBAM_MT:
    '''
    Download the necessary MT BAM and BAI from the storepath
    '''
    input:
        normal_bam = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'MT', '.bam')['store_dir'],
        normal_bai = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'MT', '.bai')['store_dir'],
        tumor_bam = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'MT', '.bam')['store_dir'],
        tumor_bai = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'MT', '.bai')['store_dir'],
    output:
        temp(alignpath_MTbqsr + '{pair}/')
    resources:
        disk = 1
    benchmark:
        benchmarkpath + 'benchmark_retrieve_BAM_MT_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Dowloading MT BAM/BAI from storepath'
    shell:
        'rsync -a {input.normal_bam} {input.normal_bai} {input.tumor_bam} {input.tumor_bai} {output}'



##############
###  MUTECT ##
##############


rule muTect_genome:
    '''
    This step uses muTect to call variants.
    '''
    input:
        rules.retrieveBAM_nuclear.output,
        muTect = muTect,
        cosmic = cosmic,
        cosmic_idx = cosmic + '.idx',
        bed = lambda wildcards: get_bed_pair(wildcards.pair),
    output:
        out = mutectpath_genome + '{pair}.tsv',
        cov = mutectpath_genome + '{pair}_cov_wig.log'
    params:
        dbsnp = dbsnp,
        ref = hg,
        normal = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'nuclear', '.bam')['wes_dir'],
        tumor = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'nuclear', '.bam')['wes_dir'],
    priority: 2
    log:
        m = mutectpath_genome_logs + '{pair}_mutect.log',
        err = mutectpath_genome_logs + '{pair}_mutect_err.log'
    conda:
        condapath + 'wes_config_conda_muTect.yaml'
    benchmark:
        benchmarkpath + 'benchmark_muTect_ref_genome_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : MuTect on genomic alignment'
    shell:
        'java -Xmx2g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --dbsnp {params.dbsnp} --intervals {input.bed} --input_file:normal {params.normal} --input_file:tumor {params.tumor} --out {output.out} --coverage_file {output.cov} > {log.m} 2> {log.err}'


rule muTect_MT:
    '''
    This step uses muTect to call MT variants.
    '''
    input:
        rules.retrieveBAM_MT.output,
        muTect = muTect,
        cosmic = cosmic,
        cosmic_idx = cosmic + '.idx',
        bed = MT_bed + '.bed',
    output:
        out = mutectpath_MT + '{pair}.tsv',
        cov = mutectpath_MT + '{pair}_cov_wig.log'
    params:
        dbsnp = dbsnp,
        ref = MT,
        normal = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'MT', '.bam')['wes_dir'],
        tumor = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'MT', '.bam')['wes_dir'],
    priority: 2
    log:
        m = mutectpath_MT_logs + '{pair}_mutect.log',
        err = mutectpath_MT_logs + '{pair}_mutect_err.log'
    conda:
        condapath + 'wes_config_conda_muTect.yaml'
    benchmark:
        benchmarkpath + 'benchmark_muTect_ref_MT_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : MuTect on MT alignment'
    shell:
        'java -Xmx2g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --dbsnp {params.dbsnp} --intervals {input.bed} --input_file:normal {params.normal} --input_file:tumor {params.tumor} --out {output.out} --coverage_file {output.cov} > {log.m} 2> {log.err}'


rule compress_mutect_cov_wig_nuclear:
    input:
        mutectpath_genome + '{pair}_cov_wig.log',
    output:
        mutectpath_genome + '{pair}_cov_wig.log.gz',
    resources:
        disk = 1
    benchmark:
        benchmarkpath + 'benchmark_compress_covwig_nuclear_ref_null_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : compressing Mutect coverage wig nuclear'
    shell:
        'gzip {input}'


rule compress_mutect_cov_wig_MT:
    input:
        mutectpath_MT + '{pair}_cov_wig.log',
    output:
        mutectpath_MT + '{pair}_cov_wig.log.gz',
    resources:
        disk = 1
    benchmark:
        benchmarkpath + 'benchmark_compress_covwig_mt_ref_null_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : compressing Mutect coverage wig MT'
    shell:
        'gzip {input}'




###############
###  VARSCAN ##
###############


rule mpileup_genome:
    '''
    This step creates an MPILEUP from a BAM.
    '''
    input:
        rules.retrieveBAM_nuclear.output,
        bed = lambda wildcards: get_bed_pair(wildcards.pair),
    output:
        temp(mpileup_genome + '{pair}/')
    params:
        ref = hg,
        normal_in = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'nuclear', '.bam')['wes_dir'],
        tumor_in = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'nuclear', '.bam')['wes_dir'],
        normal_out = lambda wildcards: get_mpileup(wildcards.pair, 'normal', 'nuclear'),
        tumor_out = lambda wildcards: get_mpileup(wildcards.pair, 'tumor', 'nuclear'),
    priority: 2
    log:
        normal = varscan_genome_logs + '{pair}_normal_mpileup.log',
        tumor = varscan_genome_logs + '{pair}_tumor_mpileup.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_mpileup_ref_genome_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Mpileup of genomic alignment'
    shell:
       'samtools mpileup -B -q 1 -f {params.ref} --positions {input.bed} {params.normal_in} > {params.normal_out} 2> {log.normal} && '
       'samtools mpileup -B -q 1 -f {params.ref} --positions {input.bed} {params.tumor_in} > {params.tumor_out} 2> {log.tumor}'


rule mpileup_MT:
    '''
    This step uses Varscan to call variants.
    '''
    input:
        rules.retrieveBAM_MT.output,
        bed = MT_bed + '.bed',
    output:
        temp(mpileup_MT + '{pair}/')
    params:
        ref = MT,
        normal_in = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'normal', 'MT', '.bam')['wes_dir'],
        tumor_in = lambda wildcards: get_bam_bai_pair(wildcards.pair, 'tumor', 'MT', '.bam')['wes_dir'],
        normal_out = lambda wildcards: get_mpileup(wildcards.pair, 'normal', 'MT'),
        tumor_out = lambda wildcards: get_mpileup(wildcards.pair, 'tumor', 'MT'),
    priority: 2
    log:
        normal = varscan_MT_logs + '{pair}_normal_mpileup.log',
        tumor = varscan_MT_logs + '{pair}_tumor_mpileup.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_mpileup_ref_MT_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Mpileup of MT alignment'
    shell:
       'samtools mpileup -B -q 1 -f {params.ref} --positions {input.bed} {params.normal_in} > {params.normal_out} 2> {log.normal} && '
       'samtools mpileup -B -q 1 -f {params.ref} --positions {input.bed} {params.tumor_in} > {params.tumor_out} 2> {log.tumor}'



rule varscan_genome:
    '''
    This step uses Varscan to call nuclear variants.
    '''
    input:
        rules.mpileup_genome.output,
        varscan = varscan,
    output:
        snv = varscanpath_genome + '{pair}_snv.tsv',
        indel = varscanpath_genome + '{pair}_indel.tsv',
    params:
        script = scriptpath + 'runVarscan.py',
        normal = lambda wildcards: get_mpileup(wildcards.pair, 'normal', 'nuclear'),
        tumor = lambda wildcards: get_mpileup(wildcards.pair, 'tumor', 'nuclear'),
        min_avg_qual = min_avg_qual,
        strand_filter = strand_filter,
        min_var_freq = min_var_freq,
        somatic_p_value = somatic_p_value,
        dataset_file = dataset_file,
        log = varscan_genome_logs + '{pair}_varscan.log'
    priority: 2
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_varscan_ref_genome_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Varscan on genomic alignment'
    script:
        '{params.script}'



rule varscan_MT:
    '''
    This step uses Varscan to call MT variants.
    '''
    input:
        rules.mpileup_MT.output,
        varscan = varscan,
    output:
        snv = varscanpath_MT + '{pair}_snv.tsv',
        indel = varscanpath_MT + '{pair}_indel.tsv',
    params:
        script = scriptpath + 'runVarscan.py',
        normal = lambda wildcards: get_mpileup(wildcards.pair, 'normal', 'MT'),
        tumor = lambda wildcards: get_mpileup(wildcards.pair, 'tumor', 'MT'),
        min_avg_qual = min_avg_qual,
        strand_filter = strand_filter,
        min_var_freq = min_var_freq,
        somatic_p_value = somatic_p_value,
        dataset_file = dataset_file,
        log = varscan_MT_logs + '{pair}_varscan.log',
    priority: 2
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_varscan_ref_MT_subject_{pair}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.pair} : Varscan on MT alignment'
    script:
        '{params.script}'



########################################
###  FINAL ANNOTATION AND FILTERING  ###
########################################

rule mutect_preAnnovar:
    input:
        genome_infile = mutectpath_genome + '{pair}.tsv',
        mt_infile = mutectpath_MT + '{pair}.tsv',
    output:
        outfile = variants_filter + '{pair}_mutect.tsv',
    params:
        script = scriptpath + 'mutectPreAnnovar.py',
        pair = '{pair}',
    message: '>> {wildcards.pair} : pre-process of Mutect for ANNOVAR'
    script:
        '{params.script}'


rule varscan_somaticFilter_preAnnovar:
    input:
        genome_snv_infile = varscanpath_genome + '{pair}_snv.tsv',
        genome_indel_infile = varscanpath_genome + '{pair}_indel.tsv',
        mt_snv_infile = varscanpath_MT + '{pair}_snv.tsv',
        mt_indel_infile = varscanpath_MT + '{pair}_indel.tsv',
    output:
        outfile = variants_filter + '{pair}_varscan.tsv',
    params:
        script = scriptpath + 'varscanSomaticFilterPreAnnovar.py',
        varscan = varscan,
        varfiltpath = variants_filter,
        pair = '{pair}',
        min_coverage = min_coverage,
        min_reads2 = min_reads2,
        min_var_freq_sf = min_var_freq_sf,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: '>> {wildcards.pair} : somaticFilter and pre-process of Varscan for ANNOVAR'
    script:
        '{params.script}'


rule merge_annotate_label:
    input:
        annovar_dbs = annovar_dbs,
        mutect_infile = variants_filter + '{pair}_mutect.tsv',
        varscan_infile = variants_filter + '{pair}_varscan.tsv',
    output:
        labelled = variants_filter + '{pair}_labelled_variants.tsv',
    params:
        script = scriptpath + 'mergeAnnotateLabel_02.py',
        tableannovar = tableannovar,
        varfiltpath = variants_filter,
        pair = '{pair}',
        humandb = humandb,
        buildver = buildver,
        buildver_mt = buildver_mt,
        protocol_genome = protocol_genome,
        protocol_mt = protocol_mt,
        operations_genome = operations_genome,
        operations_mt = operations_mt,
        rsID_dict = rsID_dict,
        processpath = processpath
#    conda:
#        condapath + 'wes_config_conda.yaml'
    message: '>> {wildcards.pair} : merging, annotating and labelling'
    script:
        '{params.script}'



###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

########################
#    GET SEVERAL FILES #
########################


rule download_dbsnp:
    '''download dbsnp_138.b37 from 1000genome '''
    output:
        dbsnp_zipped= temp(dbsnp+'.gz'),
    message: 'downloading dbsnp_138.b37 from 1000genome'
    shell:
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz && '
        'mv dbsnp_138.b37.vcf.gz {output.dbsnp_zipped}'

rule gunzip_dbsnp:
    input:
        dbsnp_zipped= dbsnp+'.gz'
    output:
        dbsnp
    message: 'unzipping dbsnp_138.b37'
    shell:
        'gunzip {input.dbsnp_zipped} || true'


# rule index_dbsnp:
#     input:
#         hg = hg,
#         dbsnp = dbsnp,
#     output:
#         dbsnp + '.idx'
#     params:
#         gatk = gatk,
#     log:
#         currentpath + '/wes_analyses/logs/index_dbsnp.log',
#     message: 'Performing ValidateVariants on dbsnp to index it'
#     conda:
#         condapath + 'wes_config_conda.yaml'
#     shell:
#         'java -jar {params.gatk} -T ValidateVariants -R {input.hg} -V {input.dbsnp} --validationTypeToExclude ALL 2> {log} || true'
#


# https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2

rule download_cosmic:
    '''download  cosmic from broadinstitute'''
    output:
        cosmic=cosmic,
    message: 'downloading b37_cosmic_v54_120711 from broadinstitute'
    shell:
        'wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf && '
        'mv b37_cosmic_v54_120711.vcf {output.cosmic}'

rule index_cosmic:
    input:
        hg = hg,
        hg_dict = hg.replace('fasta', 'dict'),
        hg_indexes = hg_indexes,
        hg_fai = hg+'.fai',
        cosmic = cosmic,
    output:
        cosmic + '.idx'
    params:
        gatk = gatk,
    log:
        processpath + 'logs/index_cosmic.log',
    message: 'Performing ValidateVariants on cosmic to index it'
    conda:
        condapath + 'wes_config_conda.yaml'
    shell:
        'java -jar {params.gatk} -T ValidateVariants -R {input.hg} -V {input.cosmic} --validationTypeToExclude ALL 2> {log} || true'



rule download_nextexp_bed:
    '''download target from illumina'''
    output:
        nextexp_bed = temp(nextexp_bed + '.bed'),
    message: 'downloading nexterarapidcapture_expandedexome_targetedregions from illumina'
    shell:
        'wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed && '
        'mv nexterarapidcapture_expandedexome_targetedregions.bed {output.nextexp_bed}'

rule fix_nextexp_bed:
    '''fix target: remove prefix chr in first column'''
    input:
        nextexp_bed + '.bed',
    output:
        nextexp_bed + '_fixed.bed',
    params:
        script = scriptpath + 'editBEDChromField_02.py',
        line_to_skip = None,
    message: 'editing nexterarapidcapture_expandedexome_targetedregions'
    script:
        '{params.script}'


rule download_nextera_bed:
    '''download target from illumina'''
    output:
        nextera_bed = temp(nextera_bed + '.bed'),
    message: 'downloading nexterarapidcapture_exome_targetedregions_v1.2 from illumina'
    shell:
        'wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed && '
        'mv nexterarapidcapture_exome_targetedregions_v1.2.bed {output}'


rule fix_nextera_bed:
    '''fix target: remove prefix chr in first column and last column(name)'''
    input:
        nextera_bed + '.bed',
    output:
        nextera_bed + '_fixed.bed',
    params:
        script = scriptpath + 'editBEDChromField_02.py',
#        line_to_skip = '0',
        line_to_skip = None,
    message: 'editing nexterarapidcapture_exome_targetedregions_v1.2'
    script:
        '{params.script}'



rule download_truseq_bed:
    '''download target from illumina'''
    output:
        truseq_bed = temp(truseq_bed + '.bed'),
    message: 'downloading truseq-exome-targeted-regions-manifest-v1-2 from illumina'
    shell:
        'wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/truseq/truseq-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'unzip truseq-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'rm truseq-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'mv truseq-exome-targeted-regions-manifest-v1-2.bed {output}'

rule fix_truseq_bed:
    '''fix target: remove prefix chr in first column and last column(name)'''
    input:
        truseq_bed + '.bed',
    output:
        truseq_bed + '_fixed.bed',
    params:
        script = scriptpath + 'editBEDChromField_02.py',
        line_to_skip = None,
    message: 'editing truseq-exome-targeted-regions-manifest-v1-2'
    script:
        '{params.script}'


rule download_truseq_rapid_bed:
    '''download target from illumina'''
    output:
        truseq_rapid_bed = temp(truseq_rapid_bed + '.bed'),
    message: 'downloading truseq-rapid-exome-targeted-regions-manifest-v1-2 from illumina'
    shell:
        'wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/truseq/truseq-rapid-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'unzip truseq-rapid-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'rm truseq-rapid-exome-targeted-regions-manifest-v1-2-bed.zip && '
        'mv truseq-rapid-exome-targeted-regions-manifest-v1-2.bed {output}'

rule fix_truseq_rapid_bed:
    '''fix target: remove prefix chr in first column and last column(name)'''
    input:
        truseq_rapid_bed + '.bed',
    output:
        truseq_rapid_bed + '_fixed.bed',
    params:
        script = scriptpath + 'editBEDChromField_02.py',
        line_to_skip = None,
    message: 'editing truseq-rapid-exome-targeted-regions-manifest-v1-2'
    script:
        '{params.script}'



rule download_nextera_dna_exome_bed:
    '''download target from illumina'''
    output:
        nextera_dna_exome_bed = temp(nextera_dna_exome_bed + '.bed'),
    message: 'downloading nextera-dna-exome-targeted-regions-manifest-v1-2 from illumina'
    shell:
        'wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/nextera-dna-exome/nextera-dna-exome-targeted-regions-manifest-bed.zip && '
        'unzip nextera-dna-exome-targeted-regions-manifest-bed.zip && '
        'rm nextera-dna-exome-targeted-regions-manifest-bed.zip && '
        'mv nextera-dna-exome-targeted-regions-manifest-v1-2.bed {output}'

rule fix_nextera_dna_exome_bed:
    '''fix target: remove prefix chr in first column and last column(name)'''
    input:
        nextera_dna_exome_bed + '.bed',
    output:
        nextera_dna_exome_bed + '_fixed.bed',
    params:
        script = scriptpath + 'editBEDChromField_02.py',
        line_to_skip = None,
    message: 'editing nextera-dna-exome-targeted-regions-manifest-v1-2'
    script:
        '{params.script}'



rule create_MTbed:
    output:
        MT_bed = MT_bed + '.bed',
    message: 'generating MT bed'
    shell:
        'echo "MT 1 16569" > {output}'

##########################################
##      DOWNLOAD ANNOVAR DATABASES      ##
##########################################

rule download_hg19refGene:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[0],
    message: 'downloading ANNOVAR hg19 refGene'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar refGene {params.humandb}'


rule download_hg19cytoBand:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[1],
    message: 'downloading ANNOVAR hg19 cytoBand'
    shell:
        '{params.annotate} -buildver hg19 -downdb cytoBand {params.humandb}'

rule download_hg19gSd:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[2],
    message: 'downloading ANNOVAR hg19 genomicSuperDups'
    shell:
        '{params.annotate} -buildver hg19 -downdb genomicSuperDups {params.humandb}'

rule download_hg19_esp:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[3],
    message: 'downloading ANNOVAR hg19 esp6500siv2_all'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar esp6500siv2_all {params.humandb}'

rule download_hg19snp138:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[4],
    message: 'downloading ANNOVAR hg19 snp138'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar snp138 {params.humandb}'

rule download_1000g2014oct:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[5],
    message: 'downloading ANNOVAR hg19 1000g2014oct_all'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar 1000g2014oct {params.humandb}'

rule download_exac03nontcga:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[6],
    message: 'downloading ANNOVAR hg19 exac03nontcga'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar exac03nontcga {params.humandb}'

rule download_ljb26_all:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[7],
    message: 'downloading ANNOVAR hg19 ljb26_all'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar ljb26_all {params.humandb}'

rule download_clinvar_20180603:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[8],
    message: 'downloading ANNOVAR hg19 clinvar_20180603'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar clinvar_20180603 {params.humandb}'

rule download_cosmic70:
    input:
        annovar,
    params:
        annotate = annotate,
        humandb = humandb,
    output:
        annovar_dbs[9],
    message: 'downloading ANNOVAR hg19 cosmic70'
    shell:
        '{params.annotate} -buildver hg19 -downdb -webfrom annovar cosmic70 {params.humandb}'



#########################################################
#                   CHECKING REQUIREMENTS               #
#########################################################


rule check_muTect:
    '''
        This step check if muTect is present in the directory set in config file.
    '''
    output:
        muTect,
    priority: 4
    shell:
        'echo "Error. muTect not found in softwares directory." && '
        'exit 1'

rule check_Annovar:
    '''
        This step check if annovar is present in the directory set in config file.
    '''
    output:
        annovar,
    priority: 3
    shell:
        'echo "Error. Annovar not found in softwares directory." && '
        'exit 1'

rule rsID_dict:
    output:
        rsID_dict,
    priority: 2
    shell:
        'echo "Error. dictionary of polymorphic DBSNP rsIDs not found." && '
        'exit 1'


############################################################
#           DEFAULT TOOLS VERSIONS NOT ON CONDA            #
############################################################

rule download_varscan2_3_9:
    output:
        varscan,
    shell:
        'wget https://downloads.sourceforge.net/project/varscan/VarScan.v2.3.9.jar && '
        'mv VarScan.v2.3.9.jar {output}'




#########################################################
#                        THE END                        #
#########################################################
