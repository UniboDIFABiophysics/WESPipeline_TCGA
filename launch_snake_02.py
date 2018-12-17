import argparse
import pandas as pd
import subprocess as sp
import yaml
import os

# get home path
homepath = os.path.expanduser('~') + '/'

# get current path
currentpath = os.getcwd() + '/'



#### Parse command line arguments and load CONFIGFILE ####

parser = argparse.ArgumentParser()
parser.add_argument('-c', required=True, dest='configfile', action='store', help='path to config file')

# get arguments
args = parser.parse_args()

### get CONFIG file absolute path
configfile = args.configfile
# in case a relative path was given as input
if not os.path.commonprefix([homepath, configfile]):
    # append to current path
    configfile = currentpath + configfile

## load config file
with open(configfile) as infile:
    config = yaml.load(infile)



# build and create PROCESSPATH & STOREPATH
processpath = homepath + config['process_path'] + config['outfolder'] + '/'
os.makedirs(processpath, exist_ok=True)

storepath = config['store_path'] + config['outfolder'] + '/'
os.makedirs(storepath, exist_ok=True)


# copy input list in PROCESSPATH
cmd = 'cp %s%s %s' %(homepath, config['input_list'], processpath)
sp.run(cmd, shell=True)


########################################################################


# set --keep-going argument
if config['keepgoing']:
    keepgoing = '--keep-going'
else:
    keepgoing = ''



# in case a relative path was given
if not os.path.commonprefix([homepath, config['snakefile']]):
    # append to homepath
    snakefile = homepath + config['snakefile']

# derive pipelinepath and change CWD to that
pipelinepath = '/'.join(snakefile.split('/')[:-1])
os.chdir(pipelinepath)

# build snakemake command
cmd = ' '.join(['snakemake',
                '-s', snakefile,
                '--configfile', configfile,
                '--use-conda',
                '--cores', config['cores'],
                '--resources disk=' + config['disk'],
                keepgoing
                ])

# run snakemake
sp.run(cmd, shell=True)

# backup every results
print('\nStoring whole analysis in storepath\n')
cmd = 'rsync -a %s %s' %(processpath, storepath)
sp.run(cmd, shell=True)


print('\n****\nEND\n****\n')
