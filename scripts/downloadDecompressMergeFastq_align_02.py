def downloadDecompressMergeFastq(name, homepath, processpath, fastq_path, fastq_cs):

    import os
    import re
    import subprocess as sp

    inpath = fastq_path.split('|')
    inoutpath = processpath + '01_fastq/01_preprocess/'
    outpath = processpath + '01_fastq/02_postprocess/'


    if not os.path.exists(inoutpath):
         os.makedirs(inoutpath)

    if not os.path.exists(outpath):
         os.makedirs(outpath)

    ##### create empty fastq1
    fastq1 = '%s_R1.fastq' %name
    cmd = '> %s%s' %(outpath, fastq1)
    sp.run(cmd, shell=True)

    ##### create empty fastq2
    fastq2 = '%s_R2.fastq' %name
    cmd = '> %s%s' %(outpath, fastq2)
    sp.run(cmd, shell=True)


    # loop over inpaths
    for ip in inpath:

        ip = homepath + ip

        # list all fastq.gz files that contain fastq_cs
        gz_list = [gz for gz in os.listdir(ip) if re.search(fastq_cs, gz)]

        # loop over them
        for gz in gz_list:

            sp.run('rsync -a %s %s' %(ip + gz, inoutpath + gz), shell=True)
            sp.run('gzip -d %s' % inoutpath + gz, shell=True)


        # list all fastq for current sample
        current_sample_fastq = [f for f in os.listdir(inoutpath) if re.search(fastq_cs, f)]

        # set all possible denominations of the paired read
        denominations = ['_R', '_sequence_', '_read']

        # loop over
        for d in denominations:

            # for current sample, list all fastq that contain current denomination
            l = [f for f in current_sample_fastq if re.search(d, f)]

            # if any exist
            if l:

                # list common substrings of each pair of fastq
                common_between_read_pairs = [d.join(re.split(d + '[12]', f)) for f in current_sample_fastq]

                # uniq list
                cbrp = sorted(list(set(common_between_read_pairs)))

                break

        # list all R1 fastq (unordered), concatenating their relative path
        R1 = [inoutpath + f for f in current_sample_fastq if re.search(d + '1', f)]

        # list all R2 fastq (unordered), concatenating their relative path
        R2 = [inoutpath + f for f in current_sample_fastq if re.search(d + '2', f)]

        # create new lists
        R1b = []
        R2b = []

        # loop over common substrings (pairs of fastqs)
        for c in cbrp:

            # re-create original string
            c = (d + '[12]').join(c.split(d))

            # loop over R1 fastq
            for r1 in R1:

                # if common substring matches R1 fastq
                if re.search(c, r1):

                    # append it to new list
                    R1b.append(r1)

            # loop over R1 fastq
            for r2 in R2:

                # if common substring matches R1 fastq
                if re.search(c, r2):

                    # append it to new list
                    R2b.append(r2)


        # join lists of ordered fastqs
        R1b = ' '.join(R1b)
        R2b = ' '.join(R2b)

        # merge current batch of fastqs to growing fastq files
        cmd = 'cat %s %s%s > %s%s' %(R1b, outpath, fastq1, inoutpath, fastq1)
        sp.run(cmd, shell=True)
        cmd = 'cat %s %s%s > %s%s' %(R2b, outpath, fastq2, inoutpath, fastq2)
        sp.run(cmd, shell=True)

        # move merged files to outpath
        cmd = 'mv %s%s %s%s' %(inoutpath, fastq1, outpath, fastq1)
        sp.run(cmd, shell=True)
        cmd = 'mv %s%s %s%s' %(inoutpath, fastq2, outpath, fastq2)
        sp.run(cmd, shell=True)

        # remove pre-processed fastq
        R = ' '.join(R1 + R2)
        sp.run('rm %s' %R, shell=True)


# Snakemake call

downloadDecompressMergeFastq(snakemake.params['name'],
                            snakemake.params['homepath'],
                            snakemake.params['processpath'],
                            snakemake.params['fastq_path'],
                            snakemake.params['fastq_cs'])
