# zsf1_gse261355
zsf1 rat kidney gse261355 pmid40409667

# mouse_adenine_rna
analysis of publicly available bulk rna seq data for

1. mouse adenine-induced CKD from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261355



Download data set
```
conda install -c bioconda sra-tools
prefetch --version
vdb-config --interactive

#!/bin/bash

dir=/mnt/h/scratch/gse261355
mkdir -P $dir && cd $dir

# Path to your input file
input_file="/mnt/h/scratch/SRR_Acc_List.txt"

# Declare an array
accession_array=()

# Read file line by line into the array
while IFS= read -r line || [[ -n "$line" ]]; do
    accession_array+=("$line")
done < "$input_file"

# Example: Print all elements
for acc in "${accession_array[@]}"; do
    echo "$acc"
done

# Loop through each accession
for acc in "${accession_array[@]}"; do
    echo "Prefetching $acc..."
    prefetch "$acc"  # download .sra file
    echo "Converting $acc to FASTQ..."
    fasterq-dump $acc/$acc.sra --split-files --outdir $acc  # convert to paired fastq
    echo "Finished $acc"
done
```

count with salmon
```
# to run locally:
SCRATCH1=/mnt/h/scratch
docker run -it --rm \
--workdir $HOME \
-v /mnt/h/scratch/gse261355:$HOME/project \
-v /mnt/g/reference:$HOME/reference \
-v $HOME:$HOME \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/h/scratch" \
-v /mnt/h/scratch:$HOME/scratch \
combinelab/salmon:latest /bin/bash

# wget -P /mnt/g/reference "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
# wget -P /mnt/g/reference "https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
# wget -P /mnt/g/reference "https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz"

# salmon index \
# -t reference/Rattus_norvegicus.GRCr8.cdna.all.fa.gz \
# -i reference/salmon_index_GRCr8 \
# --threads 4

# Path to your input file
input_file="scratch/SRR_Acc_List.txt"

# Declare an array
accession_array=()

# Read file line by line into the array
while IFS= read -r line || [[ -n "$line" ]]; do
    accession_array+=("$line")
done < "$input_file"

for acc in "${accession_array[@]}"; do
  echo "Quantifying $acc..."
  salmon quant -i reference/salmon_index_GRCr8 \
             -l A \
             -r project/$acc/"${acc}.sra.fastq" \
             --seqBias \
             --gcBias \
             --validateMappings \
             -o project/$acc/salmon_quant \
             --fldMean 200 \
             --fldSD 20 \
             -p 8
done


```




