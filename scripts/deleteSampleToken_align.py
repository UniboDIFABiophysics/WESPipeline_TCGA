def deleteSampleToken(token):

    # delete token file if --keep-going argument was specified in SNAKEMAKE
    if os.path.exists(token):
            os.remove(token)

deleteSampleToken(snakemake.params['token'])
