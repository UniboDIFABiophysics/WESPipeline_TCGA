
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


# build and create PROCESSPATH & STOREPATH
processpath = homepath + config['process_path'] + config['outfolder'] + '/'

storepath = config['store_path'] + config['outfolder'] + '/'


################################################################################

# load dataset file and set sample names as indeces
data = homepath + config['dataset']
dataset = pd.read_table(data, sep = '\t',header=0, dtype='str')
dataset = dataset.set_index('sample_id')


hg = homepath + config['fasta']['genome'] # Human Genome Reference
MT = homepath + config['fasta']['MT']
nextera = homepath + config['fasta']['nextera']
nextera_expanded = homepath + config['fasta']['nextera_expanded']
truseq = homepath + config['fasta']['truseq']
truseq_rapid = homepath + config['fasta']['truseq_rapid']
nextera_dna_exome = homepath + config['fasta']['nextera_dna_exome']

hg_indexes = [hg+'.bwt', hg+'.pac', hg+'.amb', hg+'.ann', hg+'.sa'] # Indexes for hg
MT_indexes = [MT+'.bwt', MT+'.pac', MT+'.amb', MT+'.ann', MT+'.sa']
nextera_indexes = [nextera+'.bwt', nextera+'.pac', nextera+'.amb', nextera+'.ann', nextera+'.sa']
nextexp_indexes = [nextera_expanded+'.bwt', nextera_expanded+'.pac', nextera_expanded+'.amb', nextera_expanded+'.ann', nextera_expanded+'.sa']
truseq_indexes = [truseq+'.bwt', truseq+'.pac', truseq+'.amb', truseq+'.ann', truseq+'.sa']
truseq_rapid_indexes = [truseq_rapid+'.bwt', truseq_rapid+'.pac', truseq_rapid+'.amb', truseq_rapid+'.ann', truseq_rapid+'.sa']
nextera_dna_exome_indexes = [nextera_dna_exome+'.bwt', nextera_dna_exome+'.pac', nextera_dna_exome+'.amb', nextera_dna_exome+'.ann', nextera_dna_exome+'.sa']


indels_ref = homepath+ config['ref-files']['indels_ref'] # Set of known indels
dbsnp = homepath+ config['ref-files']['dbsnp'] # SNP database


# Softwares
gatk = homepath+ config['softwares']['gatk']
adapter_removal = homepath + config['softwares']['adapter_removal']
softwares = homepath + config['folders']['softwares']
picard = homepath + config['softwares']['picard']

# Sample details
platform = config['sample-details']['platform'] # Set platform for alignment

# bed files
nextexp_bed = homepath + config['bed']['nextera_expanded'] # Set target intervals for exome analysis
nextera_bed = homepath + config['bed']['nextera']
MT_bed = homepath + config['bed']['MT']
truseq_bed = homepath + config['bed']['truseq']
truseq_rapid_bed = homepath + config['bed']['truseq_rapid']
nextera_dna_exome_bed = homepath + config['bed']['nextera_dna_exome']



# Label's parameters
n_sim = config['n_sim']
cpu_type = config['cpu_type']
thrs = config['threads']
n_cpu = config['n_cpu']

# Max threads per each multithreading rule
map_thrs = config['map_thrs']
RT_thrs = config['RT_thrs']
BaseRecal_thrs = config['BaseRecal_thrs']
PrintReads_thrs = config['PrintReads_thrs']


### GET TRIMMING PARAMS
minquality = config['adapter_removal_minquality']
minlength = config['adapter_removal_minlength']
adapters = config['adapters']
if adapters:
    pcr1 = ' --pcr1 ' + config['adapter_sequence_1']
    pcr2 = ' --pcr2 ' + config['adapter_sequence_2']
else:
    pcr1 = ''
    pcr2 = ''



#### Load input sample list #####

# move sample list file to processpath
input_list = homepath + config['input_list']

# load list of samples
with open(input_list) as infile:
    samples = infile.read().splitlines()


samples = sorted(set(samples))



# Wildcard costrains necessary for search only certain names
wildcard_constraints:
    sample = '('+'|'.join(samples)+')',


def get_kit(wildcards):
    return dataset.loc[wildcards, 'kit_wes']

def get_fastq_path(wildcards):
    path = dataset.loc[wildcards,'fastq_path_wes']
    return path

def get_fastq_cs(wildcards):
    return dataset.loc[wildcards,'fastq_common_substring_wes']


def get_bed(wildcards):
    kit = get_kit(wildcards)
    return homepath + config['bed'][kit] + '_fixed.bed'

def get_ref(wildcards):
    kit = get_kit(wildcards)
    return homepath + config['fasta'][kit]

def get_ref_indexes(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    r_indexes = [r+'.bwt', r+'.pac', r+'.amb', r+'.ann', r+'.sa']
    return r_indexes

def get_ref_fai(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    return r+'.fai'

def get_ref_dict(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    return r.replace('fasta','dict')


# define directories
tmp_dir = processpath + 'tmp/'
if not os.path.exists(tmp_dir):
     os.makedirs(tmp_dir)
preprocesspath = processpath + '01_fastq/01_preprocess/'
postprocesspath = processpath + '01_fastq/02_postprocess/'
trimmedpath = processpath + '02_fastq_trimmed/'
fastqcpath = processpath + '01_fastq/03_fastqc/'
fastqcpath_trim = processpath + '02_fastq_trimmed/01_fastqc/'

alignpath_genome = processpath + '03_alignment_genome/'
alignpath_exome = processpath + '04_alignment_exome/'
alignpath_MT = processpath + '05_alignment_MT/'

alignpath_genomebqsr = alignpath_genome + '02_bqsr/'
alignpath_MTbqsr = alignpath_MT + '02_bqsr/'
alignpath_exomeunmapped = alignpath_exome + '02_unmapped_fastq/'
fastq_logs = processpath + '01_fastq/logs/'
trimmed_logs = processpath + '02_fastq_trimmed/logs/'
genome_int = alignpath_genome + '01_intermediate/'
genome_logs = alignpath_genome + 'logs/'
exome_int = alignpath_exome + '01_intermediate/'
exome_logs = alignpath_exome + 'logs/'
MT_int = alignpath_MT + '01_intermediate/'
MT_logs = alignpath_MT + 'logs/'


storepath_trim_logs = storepath + '02_fastq_trimmed/logs/'
storepath_nuclear_bam = storepath + '03_alignment_genome/02_bqsr/'
storepath_mt_bam = storepath + '05_alignment_MT/02_bqsr/'

benchmarkpath = processpath + 'benchmarks/'

# create a folder to store tokens for NOT COMPLETED samples
if config['keepgoing']:
    ongoing = processpath + 'ongoing/'
    os.makedirs(ongoing, exist_ok=True)



##########################################
#           PIPELINE BEGINNING           #
##########################################


rule all:
    '''
      PIPELINE ENDING
    '''
    input:
        expand(fastqcpath + '{sample}_R1_fastqc.html', sample=samples),
        expand(fastqcpath + '{sample}_R2_fastqc.html', sample=samples),
        expand(fastqcpath_trim + '{sample}_R1_fastqc.html', sample=samples),
        expand(fastqcpath_trim + '{sample}_R2_fastqc.html', sample=samples),
        expand(fastqcpath + '{sample}_R1_fastqc.zip', sample=samples),
        expand(fastqcpath + '{sample}_R2_fastqc.zip', sample=samples),
        expand(fastqcpath_trim + '{sample}_R1_fastqc.zip', sample=samples),
        expand(fastqcpath_trim + '{sample}_R2_fastqc.zip', sample=samples),
        expand(trimmed_logs + '{sample}_settings.log', sample=samples),
        expand(storepath_trim_logs + '{sample}_discarded.log.gz', sample=samples),
        expand(storepath_trim_logs + '{sample}_stats.log.gz', sample=samples),
        expand(storepath_trim_logs + '{sample}_singleton.log.gz', sample=samples),
        expand(storepath_trim_logs + '{sample}_singleton_stats.log.gz', sample=samples),
        expand(storepath_nuclear_bam + '{sample}.bam', sample=samples),
        expand(storepath_nuclear_bam + '{sample}.bai', sample=samples),
        expand(storepath_mt_bam + '{sample}.bam', sample=samples),
        expand(storepath_mt_bam + '{sample}.bai', sample=samples),
    message: ' THE END '
    run:
        pass


###########################################################################################



rule downloadDecompressMergeFastq:
    output:
        outpath1 = temp(postprocesspath+ '{sample}_R1.fastq'),
        outpath2 = temp(postprocesspath+ '{sample}_R2.fastq'),
    params:
        name = '{sample}',
        homepath = homepath,
        processpath = processpath,
        fastq_path = lambda wildcards: get_fastq_path(wildcards.sample),
        fastq_cs = lambda wildcards: get_fastq_cs(wildcards.sample),
        script = scriptpath + 'downloadDecompressMergeFastq_align_02.py'
    resources:
        disk = 1
    priority: 0
    benchmark:
        benchmarkpath + 'benchmark_downloadDecompressMergeFastq_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : downloading, decompressing and merging fastq '
    script:
        '{params.script}'



rule fastqc_R1:
  input:
    fastq1 = postprocesspath+'{sample}_R1.fastq',
  output:
    fastq1_fastqc_zip = fastqcpath +'{sample}_R1_fastqc.zip',
    fastq1_fastqc_html = fastqcpath +'{sample}_R1_fastqc.html',
  params:
    outpath = fastqcpath,
  log:
    fastq_logs + '{sample}_R1_fastqc.log'
  conda:
    condapath + 'wes_config_conda.yaml'
  priority: 2
  benchmark:
      benchmarkpath + 'benchmark_fastqc1_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
  message: '>> {wildcards.sample} : fastqc1'
  shell:
    'fastqc -o {params.outpath} {input.fastq1} > {log} 2>&1'

rule fastqc_R2:
  input:
    fastq2 = postprocesspath+'{sample}_R2.fastq',
  output:
    fastq2_fastqc_zip = fastqcpath +'{sample}_R2_fastqc.zip',
    fastq2_fastqc_html = fastqcpath +'{sample}_R2_fastqc.html',
  params:
    outpath = fastqcpath,
  log:
    fastq_logs + '{sample}_R2_fastqc.log'
  conda:
    condapath + 'wes_config_conda.yaml'
  priority: 2
  benchmark:
      benchmarkpath + 'benchmark_fastqc2_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
  message: '>> {wildcards.sample} : fastqc2'
  shell:
    'fastqc -o {params.outpath} {input.fastq2} > {log} 2>&1'


rule trim:
    input:
        fastq1 = postprocesspath+'{sample}_R1.fastq',
        fastq2 = postprocesspath+'{sample}_R2.fastq',
    output:
        fastq1_trimmed = temp(trimmedpath+'{sample}_R1.fastq'),
        fastq2_trimmed = temp(trimmedpath+'{sample}_R2.fastq'),
        log_discarded = trimmed_logs + '{sample}_discarded.log',
        log_stats = trimmed_logs + '{sample}_stats.log',
        log_singleton = trimmed_logs + '{sample}_singleton.log',
        log_settings = trimmed_logs + '{sample}_settings.log',
        log_singleton_stats = trimmed_logs + '{sample}_singleton_stats.log',
    params:
        pcr1 = pcr1,
        pcr2 = pcr2,
        adapter_removal=adapter_removal,
        minquality = minquality,
        minlength = minlength,
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_trim_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Trimming'
    shell:
        '{params.adapter_removal} --file1 {input.fastq1} --file2 {input.fastq2}{params.pcr1}{params.pcr2} --stats --trimns --trimqualities --minquality {params.minquality} --minlength {params.minlength} --output1 {output.fastq1_trimmed} --output2 {output.fastq2_trimmed} --discarded {output.log_discarded} --outputstats {output.log_stats} --singleton {output.log_singleton} --settings {output.log_settings} --singletonstats {output.log_singleton_stats}'


rule compress_trim_logs:
    input:
        log_discarded = trimmed_logs + '{sample}_discarded.log',
        log_stats = trimmed_logs + '{sample}_stats.log',
        log_singleton = trimmed_logs + '{sample}_singleton.log',
        log_singleton_stats = trimmed_logs + '{sample}_singleton_stats.log',
    output:
        log_discarded = temp(trimmed_logs + '{sample}_discarded.log.gz'),
        log_stats = temp(trimmed_logs + '{sample}_stats.log.gz'),
        log_singleton = temp(trimmed_logs + '{sample}_singleton.log.gz'),
        log_singleton_stats = temp(trimmed_logs + '{sample}_singleton_stats.log.gz')
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_compress_trim_logs_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Compressing trim logs'
    shell:
        'gzip {input.log_discarded} {input.log_stats} {input.log_singleton} {input.log_singleton_stats}'


rule trim_logs_backup:
    input:
        discarded = trimmed_logs + '{sample}_discarded.log.gz',
        stats = trimmed_logs + '{sample}_stats.log.gz',
        singleton = trimmed_logs + '{sample}_singleton.log.gz',
        singleton_stats = trimmed_logs + '{sample}_singleton_stats.log.gz'
    output:
        discarded = storepath_trim_logs + '{sample}_discarded.log.gz',
        stats = storepath_trim_logs + '{sample}_stats.log.gz',
        singleton = storepath_trim_logs + '{sample}_singleton.log.gz',
        singleton_stats = storepath_trim_logs + '{sample}_singleton_stats.log.gz'
    params:
        output_dir = storepath_trim_logs
    resources:
        disk = 1
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_trim_logs_backup_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Backup of trim logs in storepath'
    shell:
        'rsync -a {input} {params.output_dir}'



rule fastqc_trimmed_R1:
  input:
    fastq1_trimmed = trimmedpath+'{sample}_R1.fastq',
  output:
    fastq1tr_fastqc_zip = fastqcpath_trim+'{sample}_R1_fastqc.zip',
    fastq1tr_fastqc_html = fastqcpath_trim+'{sample}_R1_fastqc.html',
  params:
    outpath = fastqcpath_trim,
  log:
    trimmed_logs + '{sample}_R1_fastqc.log'
  conda:
    condapath + 'wes_config_conda.yaml'
  priority: 2
  benchmark:
      benchmarkpath + 'benchmark_fastqc1_trimmed_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
  message: '>> {wildcards.sample} : fastqc1 on trimmed fastq'
  shell:
    'fastqc -o {params.outpath} {input} > {log} 2>&1'

rule fastqc_trimmed_R2:
  input:
    fastq2_trimmed = trimmedpath+'{sample}_R2.fastq',
  output:
    fastq2tr_fastqc_zip = fastqcpath_trim+'{sample}_R2_fastqc.zip',
    fastq2tr_fastqc_html = fastqcpath_trim+'{sample}_R2_fastqc.html',
  params:
    outpath = fastqcpath_trim,
  log:
    trimmed_logs + '{sample}_R2_fastqc.log'
  conda:
    condapath + 'wes_config_conda.yaml'
  priority: 2
  benchmark:
      benchmarkpath + 'benchmark_fastqc2_trimmed_ref_null_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
  message: '>> {wildcards.sample} : fastqc2 on trimmed fastq'
  shell:
    'fastqc -o {params.outpath} {input} > {log} 2>&1'



##########################
#### GENOME ALIGNMENT ####
##########################

rule map_to_genome:
    '''
    This tool maps the samples to the human reference genome.
    It does append a header necessary for the GATK analysis.
    '''
    input:
        hg_indexes,
        reference = hg,
        fastq1_trimmed = trimmedpath+'{sample}_R1.fastq',
        fastq2_trimmed = trimmedpath+'{sample}_R2.fastq',
    output:
        outfile =  temp(genome_int + '{sample}_aligned.bam'),
    params:
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = '{sample}',
    log:
        genome_logs + '{sample}_alignment.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_mapping_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Aligning to reference NUCLEAR GENOME'
    threads: map_thrs
    shell:
        'bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {input.reference} {input.fastq1_trimmed} {input.fastq2_trimmed} 2> {log} | samtools view -b > {output.outfile}'

rule sorting_genome:
    '''
    This tool sorts the input SAM or BAM file by coordinate.
    '''
    input:
        bam = genome_int + '{sample}_aligned.bam',
        picard = picard,
    output:
        outdir = temp(genome_int + '{sample}_sorted.bam'),
    params:
        tmp = tmp_dir,
    log:
        genome_logs + '{sample}_sorting.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_sorting_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : sorting genome'
    shell:
        'java -jar {input.picard}SortSam.jar INPUT={input.bam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={params.tmp} 2> {log}'

rule marking_genome:
    '''
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
    '''
    input:
        sorted_bam = genome_int + '{sample}_sorted.bam',
    output:
        out = temp(genome_int+'{sample}_marked.bam'),
    params:
        tmp = tmp_dir,
        picard = picard,
    log:
        mx = genome_logs + '{sample}_metrix.log',
        mark = genome_logs + '{sample}_marking.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_marking_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : marking genome'
    shell:
        'java -jar {params.picard}MarkDuplicates.jar'
        ' INPUT={input.sorted_bam} OUTPUT={output.out} METRICS_FILE={log.mx} '
        ' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR={params.tmp} 2> {log.mark}'

rule indexing_genome:
    '''
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    '''
    input:
        marked_bam = genome_int+'{sample}_marked.bam',
    output:
        marked_bai = temp(genome_int+'{sample}_marked.bai'),
    params:
        picard = picard,
    log:
        genome_logs + '{sample}_indexing.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_indexing_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : indexing genome'
    shell:
        'java -jar {params.picard}BuildBamIndex.jar INPUT={input.marked_bam} OUTPUT={output} 2> {log}'


rule RTC_genome:
    '''
    This tool defines intervals to target for local realignment.
    '''
    input:
        hg.replace('fasta', 'dict'),
        marked_bai = genome_int+'{sample}_marked.bai',
        ref=hg+'.fai',
        indels_ref=indels_ref,
        indels_ref_idx = indels_ref + '.idx',
        gatk = gatk,
        marked_bam = genome_int+'{sample}_marked.bam',
        bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        genome_logs+'{sample}.intervals',
    params:
        ref=hg,
    log:
        genome_logs + '{sample}_RTC.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_RTC_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{RT_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, RT_thrs=RT_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : RTC genome'
    threads: RT_thrs
    shell:
        'java -jar {input.gatk} -T RealignerTargetCreator -R {params.ref} -I {input.marked_bam} -L {input.bed} -ip 50 -known {input.indels_ref} -nt {threads} -o {output} 2> {log}'


rule IndelRealigner_genome:
    '''
    This tool performs local realignment of reads around indels.
    '''
    input:
        intvs = genome_logs+'{sample}.intervals',
        bam = genome_int+'{sample}_marked.bam',
        idx = genome_int+'{sample}_marked.bai',
    output:
        r_bam = temp(genome_int + '{sample}_realigned.bam'),
        r_idx = temp(genome_int + '{sample}_realigned.bai'),
    params:
        gatk = gatk,
        ref = hg,
        indels_ref=indels_ref,
    log:
        genome_logs + '{sample}_IndelRealigner.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_IndelRealigner_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : IndelRealigner genome'
    shell:
#        'java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.intvs} -known {params.indels_ref} -ip 50 -o {output.r_bam} 2> {log}'
        'java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.intvs} -known {params.indels_ref} -o {output.r_bam} 2> {log}'


rule BaseRecal_genome:
    '''
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    '''
    input:
        r_bam = genome_int +'{sample}_realigned.bam',
        r_idx = genome_int +'{sample}_realigned.bai',
        dbsnp = dbsnp,
        #dbsnp_idx = dbsnp + '.idx',
        bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        outtable = genome_logs + '{sample}_recal_data.table',
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    log:
        genome_logs + '{sample}_recalibrating_01.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_BaseRecal_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BaseRecal_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BaseRecal_thrs=BaseRecal_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : BaseRecal genome'
    threads: BaseRecal_thrs
    shell:
        'java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable} 2> {log}'

rule PostRecalTable_genome:
    '''
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    '''
    input:
        r_bam = genome_int +'{sample}_realigned.bam',
        r_idx = genome_int +'{sample}_realigned.bai',
        dbsnp = dbsnp,
        #dbsnp_idx = dbsnp + '.idx',
        bed = lambda wildcards: get_bed(wildcards.sample),
        outtable = genome_logs + '{sample}_recal_data.table',
    output:
        outtable_post = genome_logs + '{sample}_post_recal_data.table',
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    log:
        genome_logs + '{sample}_postrecalibrating.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_PostRecalTable_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BaseRecal_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BaseRecal_thrs=BaseRecal_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : PostRecalTable genome'
    threads: BaseRecal_thrs
    shell:
        'java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -BQSR {input.outtable} -nct {threads} -o {output.outtable_post} 2> {log}'

rule AnalyzeCovariates_genome:
    '''
    This tool creates plots to visualize base recalibration results.
    '''
    input:
        outtable1 = genome_logs + '{sample}_recal_data.table',
        outtable2 = genome_logs + '{sample}_post_recal_data.table',
    output:
        plots = genome_logs + '{sample}_recalibrationPlots.pdf',
    params:
        gatk = gatk,
        ref=hg,
    log:
        genome_logs + '{sample}_analyzecovariates.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_AnalyzeCovariates_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : AnalyzeCovariates genome'
    shell:
        'java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {input.outtable1} -after {input.outtable2} -plots {output.plots} 2> {log}'

rule PrintReads_genome:
    '''
    This tool writes out sequence read data.
    '''
    input:
        r_bam = genome_int+'{sample}_realigned.bam',
        r_idx = genome_int+'{sample}_realigned.bai',
        outtable = genome_logs + '{sample}_recal_data.table',
    output:
        bam = temp(alignpath_genomebqsr + '{sample}.bam'),
        bai = temp(alignpath_genomebqsr + '{sample}.bai'),
    params:
        gatk = gatk,
        ref=hg,
    log:
        genome_logs + '{sample}_recalibrating_02.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_PrintReads_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{PrintReads_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, PrintReads_thrs=PrintReads_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : PrintReads genome'
    threads: PrintReads_thrs
    shell:
        'java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable} -nct {threads} -o {output.bam} 2> {log}'




######################################
####  GENOME ALIGNMENT STATISTICS ####
######################################

rule samtools_bam_stats:
    '''
    Compute SAMTOOLS statistics of BAM (the final one, after PrintReads)
    '''
    input:
        alignpath_genomebqsr + '{sample}.bam',
    output:
        temp(alignpath_genomebqsr + '{sample}_samtools_bam_stats.log'),
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 4
    benchmark:
        benchmarkpath + 'benchmark_samtools_bam_stats_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Computing SAMTOOLS BAM statistics'
    shell:
        'samtools stats {input} > {output}'

rule bamToBed:
    '''
    Transform BAM to bed
    '''
    input:
        alignpath_genomebqsr + '{sample}.bam',
    output:
        temp(alignpath_genomebqsr + '{sample}.bed'),
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 4
    benchmark:
        benchmarkpath + 'benchmark_bamToBed_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Transforming BAM to BED'
    shell:
        'bedtools bamtobed -i {input} > {output}'


rule intersectBed:
    '''
    Extract reads that mapped OFF-TARGET
    '''
    input:
        bam_bed = alignpath_genomebqsr + '{sample}.bed',
        kit_bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        temp(alignpath_genomebqsr + '{sample}_off-target.bed'),
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 4
    benchmark:
        benchmarkpath + 'benchmark_intersectBed_ref_genome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Extracting OFF-TARGET reads'
    shell:
        'bedtools intersect -a {input.bam_bed} -b {input.kit_bed} -v > {output}'


rule overall_bam_stats:
    '''
    Compute and merge overall BAM statistics into one table
    '''
    input:
        bam_stats = alignpath_genomebqsr + '{sample}_samtools_bam_stats.log',
        off_target_bed = alignpath_genomebqsr + '{sample}_off-target.bed',
        kit_bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        overall_bam_stats = genome_logs + '{sample}_overall_bam_stats.log',
    params:
        script = scriptpath + 'computeOverallBAMstats.py'
    priority: 4
    message: '>> {wildcards.sample} : Computing and merging overall BAM statistics'
    script:
        '{params.script}'







##########################
####  EXOME ALIGNMENT ####
##########################

rule map_to_exome:
    '''
    This tool maps the samples to the reference exome.
    It does append a header necessary for the GATK analysis.
    '''
    input:
        hg = hg,
        hg_indexes = hg_indexes,
        fai = lambda wildcards: get_ref_fai(wildcards.sample),
        fasta_dict = lambda wildcards: get_ref_dict(wildcards.sample),
        indexes = lambda wildcards: get_ref_indexes(wildcards.sample),
        reference = lambda wildcards: get_ref(wildcards.sample),
        fastq1_trimmed = trimmedpath+'{sample}_R1.fastq',
        fastq2_trimmed = trimmedpath+'{sample}_R2.fastq',
    output:
        outfile =  temp(exome_int + '{sample}_unmapped.bam'),
    params:
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = '{sample}',
    log:
        exome_logs + '{sample}_alignment.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_mapping_ref_exome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Aligning to reference EXOME'
    threads: map_thrs
    shell:
        'bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {input.reference} {input.fastq1_trimmed} {input.fastq2_trimmed} 2> {log} | samtools view -b -f 4 > {output.outfile}'

rule sorting_exome:
    '''
    This tool sorts the input SAM or BAM file by coordinate.
    '''
    input:
        unm_bam = exome_int + '{sample}_unmapped.bam',
        picard = picard,
    output:
        outdir = temp(exome_int + '{sample}_unmapped_sorted.bam'),
    params:
        tmp = tmp_dir,
    log:
        exome_logs + '{sample}_sorting.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_sorting_ref_exome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : sorting exome'
    shell:
        'java -jar {input.picard}SortSam.jar INPUT={input.unm_bam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={params.tmp} 2> {log}'



rule extractUnmapped_convert:
    input:
        exome_int + '{sample}_unmapped_sorted.bam',
    output:
        unmapped1 = temp(alignpath_exomeunmapped + '{sample}_R1.unmapped.fastq'),
        unmapped2 = temp(alignpath_exomeunmapped + '{sample}_R2.unmapped.fastq'),
        unpair = temp(alignpath_exomeunmapped + '{sample}_unpaired.unmapped.fastq'),
    params:
        picard = picard,
    log:
        exome_logs + '{sample}_unmapped_fastq.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_Unmapped_convert_ref_exome_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : extractUnmapped convert'
    shell:
#        'java -jar -Xmx4g {params.picard}SamToFastq.jar I={input} F={output.unmapped1} F2={output.unmapped2} FU={output.unpair} VALIDATION_STRINGENCY=LENIENT 2> {log}'
        'java -jar {params.picard}SamToFastq.jar I={input} F={output.unmapped1} F2={output.unmapped2} FU={output.unpair} VALIDATION_STRINGENCY=LENIENT 2> {log}'



##########################
####   MT ALIGNMENT   ####
##########################

rule map_to_MT:
    '''
    This tool maps the samples to the mitochondria genome.
    It does append a header necessary for the GATK analysis.
    '''
    input:
        MT_indexes,
        reference = MT,
        unmapped1 = alignpath_exomeunmapped + '{sample}_R1.unmapped.fastq',
        unmapped2 = alignpath_exomeunmapped + '{sample}_R2.unmapped.fastq',
        unpair = alignpath_exomeunmapped + '{sample}_unpaired.unmapped.fastq',
    output:
        outfile =  temp(MT_int + '{sample}_aligned.bam'),
    params:
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = '{sample}',
    log:
        MT_logs + '{sample}_alignment.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 2
    benchmark:
        benchmarkpath + 'benchmark_mapping_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Aligning to reference MITOCHONDRIAL GENOME'
    threads: map_thrs
    shell:
        'bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {input.reference} {input.unmapped1} {input.unmapped2} 2> {log} | samtools view -b > {output.outfile}'

rule sorting_MT:
    '''
    This tool sorts the input SAM or BAM file by coordinate.
    '''
    input:
        bam = MT_int + '{sample}_aligned.bam',
        picard = picard,
    output:
        outdir = temp(MT_int + '{sample}_sorted.bam'),
    params:
        tmp = tmp_dir,
    log:
        MT_logs + '{sample}_sorting.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_sorting_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : sorting MT'
    shell:
        'java -jar {input.picard}SortSam.jar INPUT={input.bam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={params.tmp} 2> {log}'

rule marking_MT:
    '''
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
    '''
    input:
        sorted_bam = MT_int + '{sample}_sorted.bam',
    output:
        out = temp(MT_int+'{sample}_marked.bam'),
    params:
        tmp = tmp_dir,
        picard = picard,
    log:
        mx = MT_logs + '{sample}_metrix.log',
        mark = MT_logs + '{sample}_marking.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_marking_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : marking MT'
    shell:
        'java -jar {params.picard}MarkDuplicates.jar'
        ' INPUT={input.sorted_bam} OUTPUT={output.out} METRICS_FILE={log.mx} '
        ' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR={params.tmp} 2> {log.mark}'

rule indexing_MT:
    '''
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    '''
    input:
        marked_bam = MT_int+'{sample}_marked.bam',
    output:
        marked_bai = temp(MT_int+'{sample}_marked.bai'),
    params:
        picard = picard,
    log:
        MT_logs + '{sample}_indexing.log',
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 1
    benchmark:
        benchmarkpath + 'benchmark_indexing_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : indexing MT'
    shell:
        'java -jar {params.picard}BuildBamIndex.jar INPUT={input.marked_bam} OUTPUT={output} 2> {log}'


rule RTC_MT:
    '''
    This tool defines intervals to target for local realignment.
    '''
    input:
        MT.replace('fasta', 'dict'),
        marked_bai = MT_int+'{sample}_marked.bai',
        ref = MT+'.fai',
        indels_ref=indels_ref,
        indels_ref_idx = indels_ref + '.idx',
        gatk = gatk,
        marked_bam = MT_int+'{sample}_marked.bam',
        bed = MT_bed + '.bed',
    output:
        MT_logs+'{sample}.intervals',
    params:
        reference=MT,
    log:
        MT_logs + '{sample}_RTC.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_RTC_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{RT_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, RT_thrs=RT_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : RTC MT'
    threads: RT_thrs
    shell:
        'java -jar {input.gatk} -T RealignerTargetCreator -R {params.reference} -I {input.marked_bam} -L {input.bed} -ip 50 -known {input.indels_ref} -nt {threads} -o {output} 2> {log}'


rule IndelRealigner_MT:
    '''
    This tool performs local realignment of reads around indels.
    '''
    input:
        intvs = MT_logs+'{sample}.intervals',
        bam = MT_int+'{sample}_marked.bam',
        idx = MT_int+'{sample}_marked.bai',
    output:
        r_bam = temp(MT_int+'{sample}_realigned.bam'),
        r_idx = temp(MT_int+'{sample}_realigned.bai'),
    params:
        gatk = gatk,
        ref = MT,
        indels_ref=indels_ref,
    log:
        MT_logs + '{sample}_IndelRealigner.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_IndelRealigner_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : IndelRealigner MT'
    shell:
        'java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.intvs} -known {params.indels_ref} -ip 50 -o {output.r_bam} 2> {log}'


rule BaseRecal_MT:
    '''
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    '''
    input:
        r_bam = MT_int+'{sample}_realigned.bam',
        r_idx = MT_int+'{sample}_realigned.bai',
        dbsnp = dbsnp,
        #dbsnp_idx = dbsnp + '.idx',
        bed = MT_bed + '.bed',
    output:
        outtable = MT_logs + '{sample}_recal_data.table',
    params:
        gatk = gatk,
        ref=MT,
        indels_ref=indels_ref
    log:
        MT_logs + '{sample}_recalibrating_01.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_BaseRecal_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BaseRecal_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BaseRecal_thrs=BaseRecal_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : BaseRecal MT'
    threads: BaseRecal_thrs
    shell:
        'java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable} 2> {log}'

rule PostRecalTable_MT:
    '''
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    '''
    input:
        r_bam = MT_int +'{sample}_realigned.bam',
        r_idx = MT_int +'{sample}_realigned.bai',
        dbsnp = dbsnp,
        #dbsnp_idx = dbsnp + '.idx',
        bed = MT_bed + '.bed',
        outtable = MT_logs + '{sample}_recal_data.table',
    output:
        outtable_post = MT_logs + '{sample}_post_recal_data.table',
    params:
        gatk = gatk,
        ref=MT,
        indels_ref=indels_ref
    log:
        MT_logs + '{sample}_postrecalibrating.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_PostRecalTable_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BaseRecal_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BaseRecal_thrs=BaseRecal_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : PostRecalTable MT'
    threads: BaseRecal_thrs
    shell:
        'java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -BQSR {input.outtable} -nct {threads} -o {output.outtable_post} 2> {log}'

rule AnalyzeCovariates_MT:
    '''
    This tool creates plots to visualize base recalibration results.
    '''
    input:
        outtable1 = MT_logs + '{sample}_recal_data.table',
        outtable2 = MT_logs + '{sample}_post_recal_data.table',
    output:
        plots = MT_logs + '{sample}_recalibrationPlots.pdf',
    params:
        gatk = gatk,
        ref=MT,
    log:
        MT_logs + '{sample}_recalibrating_02.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_AnalyzeCovariates_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : AnalyzeCovariates MT'
    shell:
        'java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {input.outtable1} -after {input.outtable2} -plots {output.plots} 2> {log}'



rule PrintReads_MT:
    '''
    This tool writes out sequence read data.
    '''
    input:
        r_bam = MT_int+'{sample}_realigned.bam',
        r_idx = MT_int+'{sample}_realigned.bai',
        outtable = MT_logs + '{sample}_recal_data.table',
    output:
        bam = temp(alignpath_MTbqsr + '{sample}.bam'),
        bai = temp(alignpath_MTbqsr + '{sample}.bai'),
    params:
        gatk = gatk,
        ref=MT,
    log:
        MT_logs + '{sample}_recalibrating_02.log'
    conda:
        condapath + 'wes_config_conda.yaml'
    priority: 3
    benchmark:
        benchmarkpath + 'benchmark_PrintReads_ref_MT_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{PrintReads_thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, PrintReads_thrs=PrintReads_thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : PrintReads MT'
    threads: PrintReads_thrs
    shell:
        'java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable} -nct {threads} -o {output.bam} 2> {log}'


####################
###  BAM BACKUP  ###
####################

# rule bam_backup:
#     input:
#         n_bam = alignpath_genomebqsr + '{sample}.bam',
#         n_bai = alignpath_genomebqsr + '{sample}.bai',
#         mt_bam = alignpath_MTbqsr + '{sample}.bam',
#         mt_bai = alignpath_MTbqsr + '{sample}.bai',
#         overall_bam_stats = genome_logs + '{sample}_overall_bam_stats.log',
#         n_recal = genome_logs + '{sample}_recalibrationPlots.pdf',
#         mt_recal = MT_logs + '{sample}_recalibrationPlots.pdf',
#     output:
#         storepath_nuclear_bam + '{sample}.bam',
#         storepath_nuclear_bam + '{sample}.bai',
#         storepath_mt_bam + '{sample}.bam',
#         storepath_mt_bam + '{sample}.bai',
#     params:
#         n_output_dir = storepath_nuclear_bam,
#         mt_output_dir = storepath_mt_bam,
#         token = failedpath + '{sample}'
#     resources:
#         disk = 1
#     benchmark:
#         benchmarkpath + 'benchmark_bam_backup_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
#     message: '>> {wildcards.sample} : Backup of nuclear and MT BAM in storepath'
#     shell:
#         'rsync -a {input.n_bam} {input.n_bai} {params.n_output_dir} && '
#         'rsync -a {input.mt_bam} {input.mt_bai} {params.mt_output_dir} && '
#         'rm {params.token}'


rule bam_backup:
    input:
        bam = alignpath_genomebqsr + '{sample}.bam',
        bai = alignpath_genomebqsr + '{sample}.bai',
        overall_bam_stats = genome_logs + '{sample}_overall_bam_stats.log',
        n_recal = genome_logs + '{sample}_recalibrationPlots.pdf',
    output:
        bam = storepath_nuclear_bam + '{sample}.bam',
        bai = storepath_nuclear_bam + '{sample}.bai',
    params:
        output_dir = storepath_nuclear_bam
    resources:
        disk = 1
    priority: 4
    benchmark:
        benchmarkpath + 'benchmark_bam_backup_mt_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Backup of nuclear BAM in storepath'
    shell:
        'rsync -a {input.bam} {input.bai} {params.output_dir}'


rule bam_backup_mt:
    input:
        bam = alignpath_MTbqsr + '{sample}.bam',
        bai = alignpath_MTbqsr + '{sample}.bai',
        mt_recal = MT_logs + '{sample}_recalibrationPlots.pdf',
    output:
        bam = storepath_mt_bam + '{sample}.bam',
        bai = storepath_mt_bam + '{sample}.bai',
    params:
        output_dir = storepath_mt_bam
    resources:
        disk = 1
    priority: 4
    benchmark:
        benchmarkpath + 'benchmark_bam_backup_mt_subject_{sample}' + '_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: '>> {wildcards.sample} : Backup of MT BAM in storepath'
    shell:
        'rsync -a {input.bam} {input.bai} {params.output_dir}'


rule delete_sample_token:
    input:
        storepath_nuclear_bam + '{sample}.bam',
        storepath_nuclear_bam + '{sample}.bai',
        storepath_mt_bam + '{sample}.bam',
        storepath_mt_bam + '{sample}.bai',
#     params:
#         script = scriptpath + 'deleteSampleToken_align.py',
#         token = ongoing + '{sample}'
    message: '>> {wildcards.pair} : deleting sample token'
    run:
        token = processpath + 'ongoing/{sample}'
        if os.path.exists(token):
            os.remove(token)
    # script:
    #     '{params.script}'


###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

#################################################
#                   GET FASTA                   #
#################################################

rule download_reference:
    '''download the hg19 human reference genome from 1000genome'''
    output:
        zipped = temp(hg+'.gz'),
    version: 0.1
    benchmark:
        benchmarkpath + 'benchmark_downloadreference_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'downloading HG19'
    shell:
        'wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && '
        'mv human_g1k_v37.fasta.gz {output.zipped} '

rule gunzip_reference:
    input:
        zipped = hg+'.gz'
    output:
        hg
    benchmark:
        benchmarkpath + 'benchmark_gunzip_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'unzipping HG19'
    shell:
        'gunzip {input.zipped} || true'

rule get_nextera_fasta:
    input:
        hg = hg,
        bed_fixed = nextera_bed + '_fixed.bed'
    output:
        nextera,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating nextera fasta'
    shell:
        'bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}'

rule get_nextexp_fasta:
    input:
        hg = hg,
        bed_fixed = nextexp_bed + '_fixed.bed'
    output:
        nextera_expanded,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating nextera expanded fasta'
    shell:
        'bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}'

rule get_truseq_fasta:
    input:
        hg = hg,
        bed_fixed = truseq_bed + '_fixed.bed'
    output:
        truseq,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating truseq fasta'
    shell:
        'bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}'


rule get_truseq_rapid_fasta:
    input:
        hg = hg,
        bed_fixed = truseq_rapid_bed + '_fixed.bed'
    output:
        truseq_rapid,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating truseq_rapid fasta'
    shell:
        'bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}'


rule get_nextera_dna_exome_fasta:
    input:
        hg = hg,
        bed_fixed = nextera_dna_exome_bed + '_fixed.bed'
    output:
        nextera_dna_exome,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating nextera_dna_exome fasta'
    shell:
        'bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}'


rule get_MT_fasta:
    input:
        hg,
    output:
        MT,
    conda:
        condapath + 'wes_config_conda.yaml'
    message: 'generating MT fasta'
    shell:
        'samtools faidx {input} MT > {output}'


########################
#    GET SEVERAL FILES #
########################
rule download_indels_ref:
    '''download Mills_and_1000G_gold_standard.indels.b37 from 1000genome '''
    output:
        indel_zipped= temp(indels_ref+'.gz'),
    message: 'downloading Mills_and_1000G_gold_standard.indels.b37 from 1000genome'
    shell:
        'wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz && '
        'mv Mills_and_1000G_gold_standard.indels.b37.vcf.gz {output.indel_zipped}'

rule gunzip_indelref:
    input:
        indel_zipped = indels_ref+'.gz'
    output:
        indels_ref
    message: 'unzipping Mills_and_1000G_gold_standard.indels.b37'
    shell:
        'gunzip {input.indel_zipped} || true'

rule index_indelref:
    input:
        hg=hg,
        hg_dict = hg.replace('fasta', 'dict'),
        hg_indexes = hg_indexes,
        hg_fai = hg+'.fai',
        indel_ref = indels_ref,
    output:
        indels_ref + '.idx'
    params:
        gatk = gatk,
    log:
        processpath + 'logs/index_indels.log',
    message: 'Performing ValidateVariants on Mills_and_1000G_gold_standard.indels.b37 to index it'
    conda:
        condapath + 'wes_config_conda.yaml'
    shell:
        'java -jar {params.gatk} -T ValidateVariants -R {input.hg} -V {input.indel_ref} --validationTypeToExclude ALL 2> {log} || true'


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


#########################################################
#                   REFERENCES INDEXING                  #
#########################################################

###############
##     HG    ##
###############
rule hg_index_bwa:
    '''
    Generate the index of the reference genome for the bwa program.
    '''
    input:
        hg = hg,
        fai = hg + '.fai',
    output:
        hg_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_hg_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing hg19 with bwa'
    shell:
        'bwa index -a bwtsw {hg}'

rule index_picard_hg:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        hg = hg,
        fai = hg + '.fai',
        picard = picard,
    output:
        hg.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_hg_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing hg19 with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.hg} O={output}'

rule index_samtools_hg:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        hg = hg,
    output:
        hg + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_hg_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing hg19 with samtools'
    shell:
        'samtools faidx {input.hg} '


###############
##  Nextera  ##
###############

rule nextera_index_bwa:
    input:
        nextera = nextera,
        fai = nextera + '.fai',
    output:
        nextera_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_nextera:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        nextera = nextera,
        picard = picard,
        fai = nextera + '.fai',
    output:
        nextera.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.nextera} O={output}'

rule index_samtools_nextera:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        nextera,
    output:
        nextera + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera with samtools'
    shell:
        'samtools faidx {input} '


########################
##  Nextera_expanded  ##
########################

rule nextexp_index_bwa:
    input:
        nextera_expanded = nextera_expanded,
        fai = nextera_expanded + '.fai',
    output:
        nextexp_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_expanded_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera expanded with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_nextexp:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        nextera_expanded = nextera_expanded,
        picard = picard,
        fai = nextera_expanded + '.fai',
    output:
        nextera_expanded.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_expanded_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera expanded with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.nextera_expanded} O={output}'

rule index_samtools_nextexp:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        nextera_expanded,
    output:
        nextera_expanded + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_expanded_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera expanded with samtools'
    shell:
        'samtools faidx {input} '

##############
##  Truseq  ##
##############
rule truseq_index_bwa:
    input:
        truseq = truseq,
        fai = truseq + '.fai',
    output:
        truseq_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_truseq:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        truseq = truseq,
        picard = picard,
        fai = truseq + '.fai',
    output:
        truseq.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.truseq} O={output}'

rule index_samtools_truseq:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        truseq,
    output:
        truseq + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq with samtools'
    shell:
        'samtools faidx {input} '


###################
##  Truseq rapid ##
###################
rule truseq_rapid_index_bwa:
    input:
        truseq_rapid = truseq_rapid,
        fai = truseq_rapid + '.fai',
    output:
        truseq_rapid_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_rapid_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq_rapid with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_truseq_rapid:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        truseq_rapid = truseq_rapid,
        picard = picard,
        fai = truseq_rapid + '.fai',
    output:
        truseq_rapid.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_rapid_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq_rapid with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.truseq_rapid} O={output}'

rule index_samtools_truseq_rapid:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        truseq_rapid,
    output:
        truseq_rapid + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_truseq_rapid_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing truseq_rapid with samtools'
    shell:
        'samtools faidx {input} '


########################
##  Nextera DNA exome ##
########################
rule nextera_dna_exome_index_bwa:
    input:
        nextera_dna_exome = nextera_dna_exome,
        fai = nextera_dna_exome + '.fai',
    output:
        nextera_dna_exome_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_dna_exome_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera_dna_exome with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_nextera_dna_exome:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        nextera_dna_exome = nextera_dna_exome,
        picard = picard,
        fai = nextera_dna_exome + '.fai',
    output:
        nextera_dna_exome.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_dna_exome_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera_dna_exome with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.nextera_dna_exome} O={output}'

rule index_samtools_nextera_dna_exome:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        nextera_dna_exome,
    output:
        nextera_dna_exome + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_nextera_dna_exome_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing nextera_dna_exome with samtools'
    shell:
        'samtools faidx {input} '


#############
#     MT    #
#############
rule MT_index_bwa:
    input:
        MT = MT,
        fai = MT + '.fai',
    output:
        MT_indexes,
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_MT_index_bwa_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing MT with bwa'
    shell:
        'bwa index -a bwtsw {input}'

rule index_picard_MT:
    '''
    Generate the index of the reference genome for the picard program.
    '''
    input:
        MT = MT,
        picard = picard,
        fai = MT + '.fai',
    output:
        MT.replace('fasta', 'dict'),
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_MT_index_picard_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing MT with picard'
    shell:
        'java -jar {input.picard}CreateSequenceDictionary.jar R={input.MT} O={output}'

rule index_samtools_MT:
    '''
    Generate the index of the reference genome for the samtools and gatk programs.
    '''
    input:
        MT,
    output:
        MT + '.fai',
    conda:
        condapath + 'wes_config_conda.yaml'
    benchmark:
        benchmarkpath + 'benchmark_MT_index_samtools_ref_null_subject_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt'.format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    message: 'indexing MT with samtools'
    shell:
        'samtools faidx {input} '


#########################################################
#                   CHECKING REQUIREMENTS               #
#########################################################

rule check_AdapterRemoval:
    '''
        This step check if AdapterRemoval is present in the directory set in config file.
    '''
    output:
        adapter_removal,
    priority: 2
    shell:
        'echo "Error. AdapterRemoval not found in softwares directory." && '
        'exit 1'

rule check_GATK:
    '''
        This step check if GATK is present in the directory set in config file.
    '''
    output:
        gatk,
    priority: 5
    shell:
        'echo "Error. Genome Analysis ToolKit not found in softwares directory." && '
        'exit 1'



############################################################
#           DEFAULT TOOLS VERSIONS NOT ON CONDA            #
############################################################

rule download_picard1_119:
    output:
        picard,
    params:
        softwares=softwares,
    shell:
        'wget https://downloads.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip && '
        'unzip picard-tools-1.119.zip && '
        'rm picard-tools-1.119.zip && '
        'mv picard-tools-1.119/ {params.softwares}'




#########################################################
#                        THE END                        #
#########################################################
