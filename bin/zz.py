import pandas as pd
import re, os, subprocess

outd = "/biotron/Application/WorkStation/limengzheng/ztemp_outd/zz_tools_images"
pydir = "/biotron/Application/WorkStation/limengzheng/atomSeqTools/BaselineCreator_v1.0"


# region_df = pd.read_csv(
#     f"{outd}/zz.region.bed",
#     sep="\t",
#     index_col=None,
#     names=["chr", "sta", "end"],
# )
# with open(f"{pydir}/docs/GRCh38.microsatellites.tsv", "r") as f, open(
#     f"{outd}/zz.microsatellites.txt", "w"
# ) as fi:
#     for line in f:
#         ls = line.strip().split("\t")
#         if re.match("chromosome", line):
#             fi.write(line)
#             continue
#         for ind, row in region_df.iterrows():
#             if (
#                 row["chr"] == ls[0]
#                 and int(row["sta"]) <= int(ls[1])
#                 and int(ls[1]) <= int(row["end"])
#             ):
#                 fi.write(line)
cmd = (
    f"less -N /biotron/Application/WorkStation/limengzheng/atomSeqTools/BaselineCreator_v1.0/docs/GRCh38.microsatellites.tsv"
    f"|{{ head -n 1 ; tail -n +2|awk -v OFS='\\t' '{{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' "
    f"|bedtools intersect -a - -b "
    f"/biotron/Rawdata/Database/bed_for_analysis/T27_AtomSeq_ColonCancer_DNA_25geneAndMSI_Region_V1.1.bed  "
    f"-wa|cut -f 1,3- }} >{outd}/zz.output "
)
print(cmd)
os.system(cmd)
# subprocess.run(cmd, shell=True)
