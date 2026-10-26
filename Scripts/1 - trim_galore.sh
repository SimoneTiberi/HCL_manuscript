# Directory where fastqc files are stored:
# base-directory:
cd ~/Desktop/Emato/
mkdir trim_Filippo

# input files
fastqDir="fastq"

# BO_LYMPHO:
trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_035_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_035_100K.R2.fastq.gz

trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_035_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_035_20K.R2.fastq.gz

trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_037_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_037_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_037_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_037_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_043_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_043_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_043_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_043_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_056_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_056_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_056_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_056_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_057_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_057_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_057_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_057_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_058_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_058_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_058_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_058_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_059_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_059_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_059_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_059_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_067_100K.R1.fastq.gz $fastqDir/BO_LYMPHO_067_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/BO_LYMPHO_067_20K.R1.fastq.gz $fastqDir/BO_LYMPHO_067_20K.R2.fastq.gz


# HD:
trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_40_100K.R1.fastq.gz $fastqDir/HD_40_100K.R2.fastq.gz

trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_40_20K.R1.fastq.gz $fastqDir/HD_40_20K.R2.fastq.gz

trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_42_100K.R1.fastq.gz $fastqDir/HD_42_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_42_20K.R1.fastq.gz $fastqDir/HD_42_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_43_100K.R1.fastq.gz $fastqDir/HD_43_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_43_20K.R1.fastq.gz $fastqDir/HD_43_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_44_100K.R1.fastq.gz $fastqDir/HD_44_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_44_20K.R1.fastq.gz $fastqDir/HD_44_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_46_100K.R1.fastq.gz $fastqDir/HD_46_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_46_20K.R1.fastq.gz $fastqDir/HD_46_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_47_100K.R1.fastq.gz $fastqDir/HD_47_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_47_20K.R1.fastq.gz $fastqDir/HD_47_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_49_100K.R1.fastq.gz $fastqDir/HD_49_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_49_20K.R1.fastq.gz $fastqDir/HD_49_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_52_100K.R1.fastq.gz $fastqDir/HD_52_100K.R2.fastq.gz



trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_52_20K.R1.fastq.gz $fastqDir/HD_52_20K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_53_100K.R1.fastq.gz $fastqDir/HD_53_100K.R2.fastq.gz


trim_galore -q 20 --phred33 \
--illumina -o trim_Filippo \
--fastqc --gzip \
--paired $fastqDir/HD_53_20K.R1.fastq.gz $fastqDir/HD_53_20K.R2.fastq.gz

####################################################################################################
# in bash:
####################################################################################################
multiqc trim_Filippo --interactive -o trim_QC
