# PCAWG Consensus call image

This repo contains a docker image which encapsulates the final consensus-calling pipeline 
for SNVs and Indels for the PCAWG project.

The docker image can be built using:
```
docker build -t consensus_call .
```

## Downloading auxiliary files

The consensus calling process requires the downloading and preprocessing of several 
large standard annotation files (the reference, 1000 genomes VCF, dbSNP, repeat masker);
this is done with the following invocations:

```
mkdir dbs
docker run -it -v "${PWD}/dbs":/dbs consensus_call download reference
docker run -it -v "${PWD}/dbs":/dbs consensus_call download annotations
docker run -it -v "${PWD}/dbs":/dbs consensus_call download cosmic
```

where the cosmic downloads are optional for SNV calling (but required for indels) and
requires registration at: https://cancer.sanger.ac.uk/cosmic/register.

## Consensus ensemble calling

The step aboves creates a directory structure under `${PWD}/dbs` which the consensus pipeline
expects to be mounted at `/dbs`.  To run the conensus pipeline, invoke the docker as, eg:

```
docker run -v "${PWD}/dbs":/dbs -v "/path/to/data/files":/data consensus_call consensus snv \
    -b /data/broad_snv_file.vcf.gz \
    -d /data/dkfz_snv_file.vcf.gz \
    -m /data/muse_snv_file.vcf.gz \
    -s /data/sanger_snv_file.vcf.gz 
```

and similarly with indels (although -m is smufin). These VCFs are expected to have gone through
the PCAWG annotation pipelines (for OxoG calling and minibam-based annotation).

## Additional filtering

TODO
