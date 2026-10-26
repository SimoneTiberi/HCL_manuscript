cd ~/Desktop/Emato/

# the reference transcriptome
fasta="Kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz"

# input files
fastqDir_Filippo="trim_Filippo"

# output files
kallistoDir_Filippo="Kallisto/Filippo"

# index name
idx="Kallisto/Homo_sapiens.GRCh38.cdna.all.idx_kallisto"

# Build index
kallisto index -i $idx $fasta -t 8

# quant (to quantify transcript abundance)
# quantify trimmed reads:

# BO_LYMPHO:
kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_035_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_035_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_035_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_035_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_035_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_035_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_037_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_037_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_037_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_037_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_037_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_037_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_043_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_043_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_043_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_043_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_043_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_043_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_056_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_056_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_056_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_056_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_056_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_056_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_057_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_057_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_057_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_057_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_057_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_057_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_058_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_058_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_058_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_058_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_058_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_058_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_059_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_059_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_059_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_059_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_059_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_059_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_067_100K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_067_100K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_067_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/BO_LYMPHO_067_20K --threads 8 \
$fastqDir_Filippo/BO_LYMPHO_067_20K.R1_val_1.fq.gz $fastqDir_Filippo/BO_LYMPHO_067_20K.R2_val_2.fq.gz


# HD:
kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_40_100K --threads 8 \
$fastqDir_Filippo/HD_40_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_40_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_40_20K --threads 8 \
$fastqDir_Filippo/HD_40_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_40_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_42_100K --threads 8 \
$fastqDir_Filippo/HD_42_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_42_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_42_20K --threads 8 \
$fastqDir_Filippo/HD_42_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_42_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_43_100K --threads 8 \
$fastqDir_Filippo/HD_43_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_43_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_43_20K --threads 8 \
$fastqDir_Filippo/HD_43_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_43_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_44_100K --threads 8 \
$fastqDir_Filippo/HD_44_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_44_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_44_20K --threads 8 \
$fastqDir_Filippo/HD_44_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_44_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_46_100K --threads 8 \
$fastqDir_Filippo/HD_46_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_46_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_46_20K --threads 8 \
$fastqDir_Filippo/HD_46_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_46_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_47_100K --threads 8 \
$fastqDir_Filippo/HD_47_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_47_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_47_20K --threads 8 \
$fastqDir_Filippo/HD_47_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_47_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_49_100K --threads 8 \
$fastqDir_Filippo/HD_49_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_49_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_49_20K --threads 8 \
$fastqDir_Filippo/HD_49_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_49_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_52_100K --threads 8 \
$fastqDir_Filippo/HD_52_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_52_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_52_20K --threads 8 \
$fastqDir_Filippo/HD_52_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_52_20K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_53_100K --threads 8 \
$fastqDir_Filippo/HD_53_100K.R1_val_1.fq.gz $fastqDir_Filippo/HD_53_100K.R2_val_2.fq.gz

kallisto quant -i  $idx -o $kallistoDir_Filippo/HD_53_20K --threads 8 \
$fastqDir_Filippo/HD_53_20K.R1_val_1.fq.gz $fastqDir_Filippo/HD_53_20K.R2_val_2.fq.gz
