# aDNApipeline

Usage:

# Single ended
adnapipe_single fastq threads Population name Individual name

Example: adnapipe_single ERR2744350.fastq.gz 8 Israel I1165
# Paired
adnapipe_paired fastq_1 fastq_2 threads Population name Individual name

Example: adnapipe_paired ERR4300339_1.fastq.gz ERR4300339_2.fastq.gz 8 German kra005

The program downloads hs37d5.fa if it does not exist in the main directory. The program will also create bwa index of hs37d5.fa if it is missing.

# Notice! Do not use spaces in any arguments!
