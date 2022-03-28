source env-rnaseq/bin/activate
python3 raw_data_to_sample_sheet.py -p /mnt/tn02_bioinfo/JACA/IPS/raw_data/fastq -r /is3/projects/PUBLIC/rnaseq_anno/org/Hs/Hs.Ensembl104.fa.gz
snakemake --use-singularity --singularity-args '-B /mnt/tn02_bioinfo/JACA/IPS -B /is3/projects/PUBLIC/rnaseq_anno/org/Hs' -np
snakemake --use-singularity --singularity-args '-B /mnt/tn02_bioinfo/JACA/IPS -B /is3/projects/PUBLIC/rnaseq_anno/org/Hs' -pj70
