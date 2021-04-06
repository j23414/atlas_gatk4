#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=NX_GATK
#SBATCH --output=slurm.%x.%J.out
#SBATCH --error=slurm.%x.%J.err
#SBATCH --mail-user=jenchang@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=isu_gif_vrsc

#set -e
#set -u

start=`date +%s`

# === Load Modules here
module load singularity

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === Get input size and module versions
#module load gcc/7.3.0-xegsmw4 nextflow/20.07.1-waivkwa
NEXTFLOW=/project/isu_gif_vrsc/bin/nextflow
# === Main Program

$NEXTFLOW run 04_GATK.nf \
  --genome "00_Raw-Data/test-data/ref/*.fasta" \
  --reads "00_Raw-Data/test-data/fastq/*_{R1,R2}.fastq.gz" \
  -resume \
  -with-singularity gatk.sif \
  -with-timeline "timeline_report.html"

# $NEXTFLOW run 01_Quality-Control.nf \
#   --reads "00_Raw-Data/*.fastq.gz" \
#   -resume
#
# $NEXTFLOW run 03_GSNAP.nf \
#   --genome "02_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz" \
#   --genome_gff "02_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz" \
#   --reads "00_Raw-Data/*{1,2}.fastq.gz" \
#   -resume
#
# $NEXTFLOW run 03_Kallisto.nf \
#   --reads "00_Raw-Data/*{1,2}.fastq.gz" \
#   --genome "02_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna.gz" \
#   -resume
#$NEXTFLOW run 03_Salmon.nf \
#  --reads "00_Raw-Data/*{1,2}.fastq.gz" \
#  --genome "02_Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna.gz" \
#  -resume
# nextflow run 03_HiSat2.nf \
#   --reads "00_Raw-Data/*{R1,R2}.fastq.gz" \
#   --genome "01_Genome/*.fna.gz" \
#   --genome_gff "01_Genome/*.gff.gz" \
#   -resume
# nextflow run 01_QC-Busco.nf \
#   --genome "00_Raw-Data/*.fna" \
#   -resume
end=`date +%s`

# === Log msgs and resource use
scontrol show job ${SLURM_JOB_ID}
echo "ran NX.slurm: " `date` "; Execution time: " $((${end}-${start})) " seconds" >> LOGGER.txt
