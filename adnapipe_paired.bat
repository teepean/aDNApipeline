@echo off
title Simple aDNA Pipeline
echo.
echo            *** Simple aDNA Pipeline for paired-end fastqs ***
echo.
SETLOCAL ENABLEDELAYEDEXPANSION
if "%1"=="" goto NOPARAM

if not exist hs37d5.fa (
echo .
echo You need to have hs37d5.fa in the main directory
echo Do you want to download it?
CHOICE /C YN /M "Y/N"
IF ERRORLEVEL == 2  GOTO END
IF ERRORLEVEL == 1 (
echo Downloading and uncompressing
bin\curl.exe --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
bin\bgzip -d -@%THREADS% hs37d5.fa.gz
)
echo .
goto END)

if not exist hs37d5.fa.fai (
echo.
echo Indexing reference
bin\samtools.exe faidx hs37d5.fa
)

if not exist hs37d5.fa.bwt (
echo .
echo You need to have hs37d5.fa bwa index in the main directory
echo Do you want to index it? Notice: this might take several hours.
CHOICE /C YN /M "Y/N"
IF ERRORLEVEL == 2  GOTO END
IF ERRORLEVEL == 1 (
echo Creating an index
cygbin\bwa.exe index hs37d5.fa
)
echo .
)

set FASTQ1=%1
set FASTQ2=%2
set THREADS=%3
set POPNAME=%4
set INDNAME=%5

IF EXIST "%INDNAME%" (
  echo .
  echo Destination directory exists. Do you want to overwrite?
  echo This will delete the whole directory.
  CHOICE /C YN /M "Y/N"
  IF ERRORLEVEL == 2  (
  echo .
  echo Quitting
  echo .
  GOTO END
  mkdir %INDNAME%
  )
  IF ERRORLEVEL == 1 (
  rmdir /S /Q %INDNAME%
  mkdir %INDNAME%
  )
) ELSE (
  mkdir %INDNAME%
)

IF EXIST "%CD%\%FASTQ1%" (
set FASTQONE=%CD%\%FASTQ1%
set FASTQTWO=%CD%\%FASTQ2%
) ELSE (
  set FASTQONE=%FASTQ1%
  set FASTQTWO=%FASTQ2%
)

echo .
echo Removing adapters
echo .

cd %INDNAME%
..\cygbin\AdapterRemoval --file1 %FASTQONE% --file2 %FASTQTWO% --basename %INDNAME% --gzip --threads 2 --qualitymax 41 --collapse --preserve5p --trimns --trimqualities --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA --minlength 30 --minquality 20 --minadapteroverlap 1
..\cygbin\cat *.collapsed.gz > %INDNAME%_L0.pe.combined.tmp.fq.gz
cd ..
bin\jre\bin\java.exe -Xmx4g -jar bin\AdapterRemovalFixPrefix-0.0.5.jar %INDNAME%\%INDNAME%_L0.pe.combined.tmp.fq.gz > %INDNAME%\%INDNAME%_L0.pe.combined.fq
bin\pigz -p 2 %INDNAME%\%INDNAME%_L0.pe.combined.fq

echo .
echo Aligning
echo .


cygbin\bwa aln -t %THREADS% hs37d5.fa %INDNAME%\%INDNAME%_L0.pe.combined.fq.gz -n 0.01 -l 1024 -k 2 > %INDNAME%\%INDNAME%.sai
cygbin\bwa samse -r "@RG\tID:ILLUMINA-%INDNAME%\tSM:%INDNAME%\tPL:illumina\tPU:ILLUMINA-%INDNAME%-PE" hs37d5.fa %INDNAME%\%INDNAME%.sai %INDNAME%\%INDNAME%_L0.pe.combined.fq.gz | bin\samtools sort -@ %THREADS% -O bam - > %INDNAME%\%INDNAME%_PE.mapped.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_PE.mapped.bam

echo .
echo Dedupping
echo .

move %INDNAME%\%INDNAME%_PE.mapped.bam %INDNAME%\%INDNAME%.bam
bin\jre\bin\java.exe -Xmx4g -jar bin\DeDup-0.12.8.jar -i %INDNAME%\%INDNAME%.bam  -o %INDNAME% -u 
mv %INDNAME%\%INDNAME%.log %INDNAME%\%INDNAME%_dedup.log
bin\samtools sort -@ %THREADS% %INDNAME%\%INDNAME%_rmdup.bam -o %INDNAME%\%INDNAME%_rmdup.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_rmdup.bam

echo .
echo Damage profiling
echo .

bin\jre\bin\java.exe -Xmx4g -jar bin\DamageProfiler-0.4.9.jar -i %INDNAME%\%INDNAME%_rmdup.bam -r hs37d5.fa -l 100 -t 15 -o . -yaxis_damageplot 0.30

echo .
echo Trimming
echo .

cygbin\bam trimBam %INDNAME%\%INDNAME%_rmdup.bam %INDNAME%\tmp.bam -L 1 -R 1 
bin\samtools sort -@ %THREADS% %INDNAME%\tmp.bam -o %INDNAME%\%INDNAME%.trimmed.bam 
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%.trimmed.bam

echo .
echo Genotyping
echo .

bin\samtools mpileup -B -q 30 -Q 30 -l v42.4.1240K.pos -f hs37d5.fa %INDNAME%\%INDNAME%.trimmed.bam | bin\pileupCaller --randomHaploid   --sampleNames %INDNAME% --samplePopName %POPNAME% -f v42.4.1240K.snp -e %INDNAME%\%INDNAME%

echo .
echo Converting to plink format
echo .

echo genotypename: %INDNAME%.geno.txt > %INDNAME%\convertf.txt
echo snpname: %INDNAME%.snp.txt >> %INDNAME%\convertf.txt
echo indivname: %INDNAME%.ind.txt >> %INDNAME%\convertf.txt
echo outputformat: PACKEDPED >> %INDNAME%\convertf.txt
echo genotypeoutname: %INDNAME%.bed >> %INDNAME%\convertf.txt
echo snpoutname: %INDNAME%.bim >> %INDNAME%\convertf.txt
echo indivoutname: %INDNAME%.fam >> %INDNAME%\convertf.txt
echo outputall: YES >> %INDNAME%\convertf.txt

cd %INDNAME%
..\cygbin\convertf -p convertf.txt
cd ..

echo .
echo Done
echo The results are at %INDNAME% subdirectory
echo .

goto END
:NOPARAM
echo  Syntax:
echo     adnapipe_paired ^<fastq 1^> ^<fastq 2^> ^<threads^> ^<Population name^> ^<Individual name^>
echo.
:END


