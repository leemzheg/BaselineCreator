# BaselineCreator
Create CNV/MSI Baseline from Healthy Individuals for AtomSeq

## Requirements
- 10 or more healthy whole blood samples for CNVbaseline, 20 or more healthy whole blood samples for MSIbaseline
- Only analyse for GRCh38/hg38
- Singularity has been installed, [singularity user guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps) 
- Downloaded AtomSeqToolsDatabase
- Downloadeded AtomSeqTools's image, command:
```
singularity build AtomSeqTools_image_v2.8.sif docker://leemzheng/atomseqtools:v2.8
```
- Need a config file, fill it out like this:
```
Hg38_Fasta_Path=/PATH/to/GRCh38/hg38.fasta
Variant_library=/PATH/to/AtomSeqToolsDatabase/Variant_library
```
- If you have established AtomSeqTools hg38 alignment index before, ignore this step; else, using the command line:
```
python3 AtomSeqTools/bin/make_index.py \
-image atomSeqTools_images_v2.8.sif \
-fasta /PATH/to/GRCh38/hg38.fasta
```

## Mode one(make CNV baseline):
```
python3 BaselineCreator.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0322-LXQ-T11-A1 0327-MHT-T11-A1 0514-YZQ-T11-A1 0419-WGT-T11-A1 0408-ZKY-T11-A1   \
-outdir /PATH/to/Output \
-config configure \
-bed docs/T11_AtomSeq_LungCancer_DNA_17gene_Primer_V1.0.bed \
-baseline-type CNV -baseline-name CNVbaseline_LungCancer_17gene.cnn
```

## Mode two(make MSI baseline):
```
python3 BaselineCreator.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0314-CJQ-T27-A1 0314-WSL-T27-A1 0314-CQL-T27-A1 0314-XJ-T27-A1 0519-CSH-T27-A1 0519-CSY-T27-A1 \
-outdir /PATH/to/Output \
-config configure \
-bed docs/T27_AtomSeq_ColonCancer_DNA_25geneAndMSI_Primer_V1.1.bed \
-baseline-type MSI -baseline-name MSIbaseline_ColonCancer_25gene.msisenser-pro
```

## Contacts
If you have any questions or feedback, please contact us at: \
**Email:** mengzheng-li@ebiotron.com
