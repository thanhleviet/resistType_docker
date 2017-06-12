# resistType_docker

1. Introduction

ResistType is a pipeline to make predictions of resistance genotypes and phenotypes for bacterial species, particularly Enterobacteriaceae, from Whole Genome Sequencing (WGS) data. It uses a combination of the assembly approach and the mapping approach to make genotype predictions, matching sequences with the reference database of resistance genes using BLASTn, reporting closest match at protein sequence level where applicable. It also estimates copy number (CN) of resistance genes given that a WGS assembly from SPADes is provided where the contig names contain the contig coverage information. CN estimates are provided as an additional information for researchers to further investigate the possible effect of CNs of certain resistance genes to the phenotypic responses. In the current version of resistType, CN values of TEM is used to infer resistance to co-amoxiclav (threshold 2.5) but this needs to be intepreted with care and the threshold definition still needs additional work to finalise. 

The resistance phenotypic prediction is performed based on a catalogue of known genotype-phenotype correspondence originally compiled by Dr. Nicole Stoesser (Stoesser, N., et al. Predicting antimicrobial susceptibilities for Escherichia coli and Klebsiella pneumoniae isolates using whole genomic sequence data. The Journal of antimicrobial chemotherapy 2013;68(10):2234-2244 ) from the literature and updated on a regular basis to reflect changes to the intepretations in the field. 

The resistType package also includes some additional features that benefit users analysing Enterobacteriacaea spp. and other bacterial species in general. These features include MLST typing (for a handful of bacterial species : Escherichia coli, Enterobacter cloacae, Klebsiella pneumoniae, Klebsiella oxytoca, Pseudomonas aeruginoas, Staphyloccocus aureus, Streptococcus suis, etc. ), plasmid replication genes (inc genes) typing using the plasmidFinder database, and insertion sequence finding using insertion sequences downloaded from the ISFinder database (https://www-is.biotoul.fr). MLST is performed blinded by the species information and the typing script will pick up the species' alleles automatically from the MLST allele database. These features could be run either from raw reads or from assembled contigs in reasonable time (e.g. for a normal single culture isolate, ~ 30'-1 min/sample from raw reads, and <30'/sample from assembled contigs).  


2. Usages

2.1  Building the docker image of resistType

To create the resistType docker image, download the Dockerfile and use

> docker build . 

2.2  Using resistType image

Once the image is created with an image id, for example c9d18d875199, resistType can be run as below 

>  docker  run --rm -it  -v /Path/To/LocalDirectory/:/data c9d18d875199 plasmidType_v0.1.py -s sampleID -c assemblies/sampleID.contigs.fasta -o outDir

The input could also be raw fastq files, or bam file 

>  docker  run --rm -it  -v /Path/To/LocalDirectory/:/data c9d18d875199 plasmidType_v0.1.py -s sampleID -1 x.fq1.gz -2 x.fq2.gz -o outDir -m 1

> docker  run --rm -it  -v /Path/To/LocalDirectory/:/data c9d18d875199 plasmidType_v0.1.py -s sampleID -b x.bam -o outDir -m 1

"-m 1" option is for quick scan of resistance genes where assemblies are not created. If '-m 1' is not used, resistType will attempt to use SPAdes to assemble the raw reads, which would be time consuming. When the input is a bam file, the fastq files are generated hence it is more time consuming. 



To run MLST typing, 

> docker  run --rm -it  -v /Path/To/LocalDirectory/:/data c9d18d875199 mlstTyping.py -s sampleID -1 x.fq1.gz -2 x.fq2.gz -o outDir 

or 
> docker  run --rm -it  -v /Path/To/LocalDirectory/:/data c9d18d875199 mlstTyping.py -s sampleID -c contig.fasta -o outDir 

Similarly, use plasmidFinder.py for plasmidFinder typing, and isfinder.py for insertion sequence identification. The input of MLST, plasmidFinder and isFinder is similar. The default output directory of these are mlst/, plasmidfinderOutput/, isfinderOutput/ respectively, and all output files have a prefix that is the sampleID. 





