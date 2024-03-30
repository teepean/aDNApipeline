@echo off
title Simple aDNA Pipeline
echo.
echo            *** Simple aDNA Pipeline for single-end fastqs***
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

set FASTQ=%1
set THREADS=%2
set POPNAME=%3
set INDNAME=%4

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

IF EXIST "%CD%\%FASTQ%" (
set FASTQ1=%CD%\%FASTQ%
) ELSE (
  set FASTQ1=%FASTQ%
)

echo Preprocessing and removing adapters
bin\fastp --in1 %FASTQ1% --out1  %INDNAME%\%INDNAME%.fastp.fastq.gz --thread 4 --length_required 25 --json %INDNAME%\%INDNAME%.fastp.json --html %INDNAME%\%INDNAME%.fastp.html

bwamsvc\bwa mem -t %THREADS% -k 19 -r 2.5 -R "@RG\tID:ILLUMINA-%INDNAME%\tSM:%INDNAME%\tPL:illumina\tPU:ILLUMINA-%INDNAME%-SE" hs37d5.fa %INDNAME%\%INDNAME%.fastp.fastq.gz | bin\samtools sort --no-PG -@ %THREADS% -m2G -O bam - > %INDNAME%\%INDNAME%_SE.mapped.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_SE.mapped.bam
bin\samtools view -b %INDNAME%\%INDNAME%_SE.mapped.bam Y > %INDNAME%\%INDNAME%_Y.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_Y.bam

echo .
echo Marking duplicates
echo .

bin\jre\bin\java.exe -Xmx4g -jar bin\picard\picard.jar MarkDuplicates INPUT=%INDNAME%\%INDNAME%_SE.mapped.bam OUTPUT=%INDNAME%\%INDNAME%_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE=%INDNAME%\%INDNAME%_rmdup.metrics VALIDATION_STRINGENCY=SILENT
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_rmdup.bam

echo .
echo Trimming
echo .

cygbin\bam trimBam %INDNAME%\%INDNAME%_rmdup.bam %INDNAME%\tmp.bam -L 1 -R 1 
bin\samtools sort --no-PG -@ %THREADS% %INDNAME%\tmp.bam -o %INDNAME%\%INDNAME%.trimmed.bam 
del %INDNAME%\tmp.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%.trimmed.bam

echo .
echo Genotyping
echo .

bin\samtools mpileup -B -q 30 -Q 30 -l v42.4.1240K.pos -f hs37d5.fa %INDNAME%\%INDNAME%.trimmed.bam | bin\pileupCaller --randomHaploid   --sampleNames %INDNAME% --samplePopName %POPNAME% -f v42.4.1240K.snp -p %INDNAME%\%INDNAME%

echo .
echo Done
echo The results are at %INDNAME% subdirectory
echo .

goto END
:NOPARAM
echo  Syntax:
echo     adnapipe_single ^<fastq1^> ^<threads^> ^<Population name^> ^<Individual name^>
echo.
:END

