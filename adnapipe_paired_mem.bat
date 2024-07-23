@echo off
title Simple aDNA Pipeline
echo.
echo            *** Simple aDNA Pipeline for paired-end fastqs***
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
goto END)

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

echo
echo Preprocessing and removing adapters
echo

cd %INDNAME%
..\cygbin\adapterremoval3.exe --file1 %FASTQONE% --file2 %FASTQTWO% --basename %INDNAME% --gzip --threads %THREADS%  --qualitymax 41 --collapse --preserve5p --trimns --trimqualities --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA --minlength 30 --minquality 20 --minadapteroverlap 1 --merge --interleaved-output
copy /b *.gz %INDNAME%_full.fastq.gz
cd ..

echo
echo "Aligning"
echo

bin\bwa mem -t %THREADS% -R "@RG\tID:ILLUMINA-%INDNAME%\tSM:%INDNAME%\tPL:illumina\tPU:ILLUMINA-%INDNAME%-PE" hs37d5.fa %INDNAME%\%INDNAME%_full.fastq.gz | bin\samtools sort --no-PG -@ %THREADS% -m2G -O bam - > %INDNAME%\%INDNAME%_PE.mapped.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_PE.mapped.bam
bin\samtools view -@ %THREADS% --no-PG -b %INDNAME%\%INDNAME%_PE.mapped.bam Y > %INDNAME%\%INDNAME%_Y.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_Y.bam

echo .
echo Marking duplicates
echo .

move %INDNAME%\%INDNAME%_PE.mapped.bam %INDNAME%\%INDNAME%.bam
bin\jre\bin\java.exe -Xmx16g -jar bin\picard\picard.jar MarkDuplicates INPUT=%INDNAME%\%INDNAME%_SE.mapped.bam OUTPUT=%INDNAME%\%INDNAME%_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE=%INDNAME%\%INDNAME%_rmdup.metrics VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=500000
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_rmdup.bam

echo .
echo Genotyping
echo .

bin\samtools mpileup -B -q 30 -Q 30 -l v42.4.1240K.pos -f hs37d5.fa %INDNAME%\%INDNAME%_rmdup.bam | bin\pileupCaller --randomHaploid   --sampleNames %INDNAME% --samplePopName %POPNAME% -f v42.4.1240K.snp -p %INDNAME%\%INDNAME%  > %INDNAME%\%INDNAME%.stats.txt 2>&1

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


