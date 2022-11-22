# all of the outputs created by the below commands are available at https://osf.io/6k74c
# sourmash version 4.4.2 was used to create the test data (I think).
# if you have conda and mamba installed on your system, you can install this version of sourmash with:

mamba env create -n sourmash sourmash=4.4.2
conda activate sourmash

###################################################################################
## DOWNLOAD THE DATA
###################################################################################

# SRR5936131
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/001/SRR5936131/SRR5936131_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/001/SRR5936131/SRR5936131_2.fastq.gz

# SRR5947006
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/006/SRR5947006/SRR5947006_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/006/SRR5947006/SRR5947006_2.fastq.gz

# SRR5935765
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/005/SRR5935765/SRR5935765_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/005/SRR5935765/SRR5935765_2.fastq.gz

# SRR5936197
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/007/SRR5936197/SRR5936197_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR593/007/SRR5936197/SRR5936197_2.fastq.gz

# SRR5946923
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/003/SRR5946923/SRR5946923_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/003/SRR5946923/SRR5946923_2.fastq.gz

# SRR5946920
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/000/SRR5946920/SRR5946920_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR594/000/SRR5946920/SRR5946920_2.fastq.gz


###################################################################################
## SKETCH THE DATA
###################################################################################

for infile in *_1.fastq.gz
do
    bn=$(basename ${infile} _1.fastq.gz)
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.sig ${infile} ${bn}_2.fastq.gz
done

###################################################################################
## RUN SOURMASH GATHER
###################################################################################

wget -O gtdb-rs207.genomic-reps.dna.k31.zip https://osf.io/3a6gn/download
for infile in *sig
do
    bn=$(basename $infile .sig)
    sourmash gather ${infile} gtdb-rs207.genomic-reps.dna.k31.zip -o ${bn}_gather_gtdbrs207_reps.csv
done

###################################################################################
## RUN SOURMASH TAXONOMY
###################################################################################

wget -O gtdb-rs207.taxonomy.csv.gz https://osf.io/v3zmg/download
gunzip gtdb-rs207.taxonomy.csv.gz
sourmash tax prepare -t gtdb-rs207.taxonomy.csv -o gtdb-rs207.taxonomy.sqldb -F sql

for infile in *_gather_gtdbrs207_reps.csv
do
    sourmash tax annotate -g ${infile} -t gtdb-rs207.taxonomy.sqldb 
done

###################################################################################
## DOWNSAMPLE THE SIGNATURES TO MAKE THEM SMALLER
###################################################################################

for infile in *sig
do
    sourmash sig downsample --scaled 200000 -o ${infile} ${infile}
done


###################################################################################
## COMPRESS CSV FILES 
###################################################################################

for infile in *csv
do
   gzip ${infile}
done

####################################################################################
## RUN SOURMASH COMPARE
####################################################################################

sourmash compare -k 31 -o gut_compare.csv *sig
