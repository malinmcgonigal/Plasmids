# Summary of Steps 
1. Trim Galore: Adapter trimming and quality control of raw reads.
2. MetaSPAdes: Assembly of metagenomic sequences.
3. PlasmidFinder: Annotation of plasmids in the contigs.

Additional steps that were later discarded: 
1. MetaBAT: Binning of contigs, though not many bins were formed.
2. Kraken2: Attempted annotation of plasmids, but the attempt was unsuccessful.

# Trim Galore # 
Used trim galore for adapter trimming and quality control. 
```
conda install -c bioconda trim-galore
nano run_trim_galore.sh

#!/bin/bash

# Define arrays for sample identifiers and paths
SAMPLES=("M0388_3months/child" "M0388_3months/mother" "M0399_3months/mother" "M0399_3months/child" "M0487_3months/child" "M0487_3months/mother")
BASE_DIR="/mnt/malin_m/metagenome"

# Step 1: Run Trim Galore for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  RAW_READS_1="$BASE_DIR/raw_reads/${SAMPLE}_1.fastq"
  RAW_READS_2="$BASE_DIR/raw_reads/${SAMPLE}_2.fastq"
  OUTPUT_DIR="$BASE_DIR/trimmed/$SAMPLE"

  mkdir -p "$OUTPUT_DIR"

  trim_galore --paired "$RAW_READS_1" "$RAW_READS_2" -o "$OUTPUT_DIR"
done

chmod +x run_trim_galore.sh
./run_trim_galore.sh
```



# METASPADES
Used metaspades to assemle the sequences. Originally tried the --metaplasmid option, but found that it was more effective to use the defaul --meta option, and then use plasmidfinder to isolate plasmids. 
```
conda install -c bioconda spades
nano run_metaspades.sh

#!/bin/bash

# Define arrays for sample identifiers and paths
SAMPLES=("M0388_3months/child" "M0388_3months/mother" "M0399_3months/mother" "M0399_3months/child" "M0487_3months/child" "M0487_3months/mother")
BASE_DIR="/mnt/malin_m/metagenome"

# Step 1: Run MetaSPAdes for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  TRIMMED_READS_1="$BASE_DIR/trimmed/$SAMPLE/SRR*_1_val_1.fq"
  TRIMMED_READS_2="$BASE_DIR/trimmed/$SAMPLE/SRR*_2_val_2.fq"
  OUTPUT_DIR="$BASE_DIR/spades/${SAMPLE}_meta"

  spades.py --meta -1 "$TRIMMED_READS_1" -2 "$TRIMMED_READS_2" -o "$OUTPUT_DIR"
done

chmod +x run_metaspades.sh
./run_metaspades.sh

```

# Plasmid Finder
Used plasmid finder to annotate the plasmids against the contigs 
```

git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder
pip install -r requirements.txt

nano run_plasmidfinder.sh

#!/bin/bash

# Define arrays for sample identifiers and paths
SAMPLES=("M0388_3months/child" "M0388_3months/mother" "M0399_3months/mother" "M0399_3months/child" "M0487_3months/child" "M0487_3months/mother")
BASE_DIR="/mnt/malin_m/metagenome"
PLASMIDFINDER_DB="/path/to/genomicepidemiology-plasmidfinder_db-3e77502f23c4"

# Step 1: Run PlasmidFinder for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  CONTIGS_PATH="$BASE_DIR/spades/$SAMPLE/contigs.fasta"
  OUTPUT_DIR="$BASE_DIR/plasmidfinder/$SAMPLE"

  mkdir -p "$OUTPUT_DIR"

  python plasmidfinder.py -i "$CONTIGS_PATH" -o "$OUTPUT_DIR" -p "$PLASMIDFINDER_DB"
done

chmod +x run_plasmidfinder.sh
./run_plasmidfinder.sh
```

# METABAT 
Attempted to bin the contigs using metabat, not many bins were formed so were not used in the end for further processing. 
```
conda install -c bioconda metabat2

nano run_metabat.sh


#!/bin/bash

# Define arrays for sample identifiers and paths
SAMPLES=("M0388_3months/child" "M0388_3months/mother" "M0399_3months/mother" "M0399_3months/child" "M0487_3months/child" "M0487_3months/mother")
BASE_DIR="/mnt/malin_m/metagenome"

# Step 1: Bowtie2-build for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  CONTIGS_PATH="$BASE_DIR/spades/$SAMPLE/contigs.fasta"
  INDEX_NAME="contigs_$SAMPLE"
  
  bowtie2-build "$CONTIGS_PATH" "$INDEX_NAME"
done

# Step 2: Bowtie2 align and SAMtools sort/index for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  TRIMMED_READS_1="$BASE_DIR/trimmed/$SAMPLE/SRR*_1_val_1.fq"
  TRIMMED_READS_2="$BASE_DIR/trimmed/$SAMPLE/SRR*_2_val_2.fq"
  INDEX_NAME="contigs_$SAMPLE"
  SAM_OUTPUT="$BASE_DIR/metabat/$SAMPLE/mapped_reads.sam"
  BAM_OUTPUT="$BASE_DIR/metabat/$SAMPLE/mapped_reads.bam"
  SORTED_BAM_OUTPUT="$BASE_DIR/metabat/$SAMPLE/mapped_reads_sorted.bam"
  
  bowtie2 -x "$INDEX_NAME" -1 "$TRIMMED_READS_1" -2 "$TRIMMED_READS_2" -S "$SAM_OUTPUT"
  
  # Convert SAM to BAM, sort and index
  samtools view -bS "$SAM_OUTPUT" > "$BAM_OUTPUT"
  samtools sort "$BAM_OUTPUT" -o "$SORTED_BAM_OUTPUT"
  samtools index "$SORTED_BAM_OUTPUT"
  
  # Cleanup intermediate files
  rm "$SAM_OUTPUT" "$BAM_OUTPUT"

  # Step 3: Run MetaBAT
  runMetaBat.sh "$BASE_DIR/spades/$SAMPLE/contigs.fasta" "$SORTED_BAM_OUTPUT"
done

chmod +x process_samples.sh
./process_samples.sh
```

# Kraken2 
Attempted to use Kraken2 for annotating the plasmids, did not work.
```

git clone https://github.com/DerrickWood/kraken2.git
cd kraken2

./install_kraken2.sh .
kraken2-build --download-taxonomy --db /path/to/kraken2_db
kraken2-build --download-library plasmid --db /path/to/kraken2_db
kraken2-build --build --db /path/to/kraken2_db

nano run_kraken.sh


#!/bin/bash

# Define arrays for sample identifiers and paths
SAMPLES=("M0388_3months/child" "M0388_3months/mother" "M0399_3months/mother" "M0399_3months/child" "M0487_3months/child" "M0487_3months/mother")
BASE_DIR="/mnt/malin_m/metagenome"
KRAKEN_DB="/mnt/kohei/tools/kraken2/kraken2_db"

# Step 1: Run Kraken2 for all samples
for SAMPLE in "${SAMPLES[@]}"; do
  CONTIGS_PATH="$BASE_DIR/spades/$SAMPLE/contigs.fasta"
  OUTPUT_DIR="$BASE_DIR/kraken/$SAMPLE"
  OUTPUT_FILE="$OUTPUT_DIR/kraken2_plasmid_output.txt"
  REPORT_FILE="$OUTPUT_DIR/kraken2_plasmid_report.txt"
  
  mkdir -p "$OUTPUT_DIR"
  
  /home/malin_m/software/kraken2/kraken2 --db "$KRAKEN_DB" --output "$OUTPUT_FILE" --report "$REPORT_FILE" "$CONTIGS_PATH"
done

chmod +x run_kraken.sh
./run_kraken.sh
```
