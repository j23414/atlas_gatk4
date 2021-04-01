# GATK4 module on Atlas HPC

Testing gatk4 module on Atlas HPC


## Example Dataset

Reuse dataset from GATK nextflow pipeline ([GitHub:HuffordLab/Maize_WGS_Build](https://github.com/HuffordLab/Maize_WGS_Build))

```
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
ls -1 test-data/
#> fastq          # <= folder of paired reads
#> read-group.txt
#> ref            # <= folder containing one genome reference
```

## Load and link GATK module

**(1) Load gatk module**

```
salloc -A isu_gif_vrsc
module load singularity
module use  /project/users/jim.coyle/modulefiles
module load gatk
```

Turns out `gatk` is a function:

```
which gatk
gatk is a function
gatk () 
{ 
    singularity exec /project/users/jim.coyle/containers/gatk/4.2.0.0//gatk-4.2.0.0.sif gatk "$@"
}
```

Bind to current working directory or it will end up with strange pathing errors (can't find input files, any output files tend to be sent to home directory instead of current directory)

```
gatk2 () 
{ 
    singularity exec --bind $PWD /project/users/jim.coyle/containers/gatk/4.2.0.0//gatk-4.2.0.0.sif gatk "$@"
}
```

**(2) Print help message**

```
gatk2 --list

#> Using GATK jar /gatk/gatk-package-4.2.0.0-local.jar
#> Running:
#>     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.2.0.0-local.jar --help
#> USAGE:  <program name> [-h]
#> 
#> Available Programs:
#> --------------------------------------------------------------------------------------
#> Base Calling:                                    Tools that process sequencing machine data, e.g. Illumina base calls, and detect sequencing level attributes, e.g. adapters
#>     CheckIlluminaDirectory (Picard)              Asserts the validity for specified Illumina basecalling data.  
#>     CollectIlluminaBasecallingMetrics (Picard)   Collects Illumina Basecalling metrics for a sequencing run.  
#>     CollectIlluminaLaneMetrics (Picard)          Collects Illumina lane metrics for the given BaseCalling analysis directory.  
#>     ExtractIlluminaBarcodes (Picard)             Tool determines the barcode for each read in an Illumina lane.  
#>     IlluminaBasecallsToFastq (Picard)            Generate FASTQ file(s) from Illumina basecall read data. 
# ...
```

## Start GATK pipeline

```
# === Input
REF=00_Raw-Data/test-data/ref/b73_chr1_150000001-151000000.fasta
R1=00_Raw-Data/test-data/fastq/BioSample01_R1.fastq.gz
R2=00_Raw-Data/test-data/fastq/BioSample01_R2.fastq.gz
READNAME=BioSample01
```

The steps are slighly modified from the GATK workbook

![](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/assets/gatk-workflow.png)

**(1) FastqToSam**

```
gatk2 FastqToSam \
  --FASTQ $R1  \
  --FASTQ2 $R2  \
  --OUTPUT ${READNAME}.bam \
  --READ_GROUP_NAME ${READNAME} \
  --SAMPLE_NAME ${READNAME}_name \
  --LIBRARY_NAME ${READNAME}_lib \
  --PLATFORM ILLUMINA \
  --SEQUENCING_CENTER ISU
```

**(2) MarkIlluminaAdapters**

```
gatk2 MarkIlluminaAdapters \
  ...
```

**(2b) bwa index**

Needed to load `bwa` or `bwa-mem2` separately, was not included in `gatk` module. (We should include it, along with `samtools`)

```
bwa-mem2 index $REF
```

**(3) SamToFastq**

```
gatk2 SamToFastq \
 ...
```

**(3b) CreateSequenceDictionary**

**(4) bwa mem**

```
bwa-mem2 mem -t 16 $REF $R1 $R2 |\
  samtools view --threads 16 -bS - > ${READNAME}_aln.bam
```

**(5) MergeBamAlignment**

```
gatk2 MergeBamAlignment \
  ...
```

**(6) MarkDuplicates**

```
gatk2 MarkDuplicates \
 ...
```

**(7) AddOr-ReplaceReadGroups**

...





