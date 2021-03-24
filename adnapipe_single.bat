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
bin\curl.exe -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
bin\bgzip -d hs37d5.fa.gz
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
  rmdir /S %INDNAME%
  mkdir %INDNAME%
  ) 
) ELSE (
  mkdir %INDNAME%
)


echo .
echo Removing adapters
echo .

IF EXIST "%CD%\%FASTQ%" (
set FASTQ1=%CD%\%FASTQ%
) ELSE (
  set FASTQ1=%FASTQ%
)

cd %INDNAME%

..\cygbin\AdapterRemoval --file1 %FASTQ1% --basename %INDNAME% --gzip --threads 2 --qualitymax 41 --collapse --preserve5p --trimns --trimqualities --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA --minlength 30 --minquality 20 --minadapteroverlap 1
cd ..
echo .
echo Aligning
echo .


cygbin\bwa aln -t %THREADS% hs37d5.fa %INDNAME%\%INDNAME%.truncated.gz -n 0.01 -l 1024 -k 2 > %INDNAME%\%INDNAME%.sai
cygbin\bwa samse -r "@RG\tID:ILLUMINA-%INDNAME%\tSM:%INDNAME%\tPL:illumina\tPU:ILLUMINA-%INDNAME%-SE" hs37d5.fa %INDNAME%\%INDNAME%.sai %INDNAME%\%INDNAME%.truncated.gz | bin\samtools sort -@ %THREADS% -O bam - > %INDNAME%\%INDNAME%_SE.mapped.bam
bin\samtools index -@ %THREADS% %INDNAME%\%INDNAME%_SE.mapped.bam

echo .
echo Marking duplicates
echo .

move %INDNAME%\%INDNAME%_SE.mapped.bam %INDNAME%\%INDNAME%.bam
bin\jre\bin\java.exe -Xmx4g -jar bin\picard\picard.jar MarkDuplicates INPUT=%INDNAME%\%INDNAME%.bam OUTPUT=%INDNAME%\%INDNAME%_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE=%INDNAME%\%INDNAME%_rmdup.metrics VALIDATION_STRINGENCY=SILENT
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
echo     adnapipe_single ^<fastq1^> ^<threads^> ^<Population name^> ^<Individual name^>
echo.
:END


