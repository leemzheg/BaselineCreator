#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/05/20 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os, glob
import argparse, textwrap
import pandas as pd
import datetime, yaml
from pathlib import Path

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="Create CNV/MSI Baseline from Healthy Individuals",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="The bed file of primers should contain at least first six columns(delimited by table):\n"
        "1. Chormosome ID  (consistent with reference genome);\n"
        "2. Start of primer (0-based);\n"
        "3. End of primer (1-based);\n"
        "4. Primer ID;\n"
        "5. Primer length;\n"
        "6. Targeted strand (+/-).",
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-image", metavar="IMAGE", help="    Singularity image [*.sif]", required=True
    )
    required.add_argument(
        "-fq-dir",
        metavar="PATH",
        help="    The path to directory of raw fastq",
        required=True,
    )
    required.add_argument(
        "-fq-prefix",
        metavar="STR",
        help="    Specify fastq prefix to be analyzed(can be many, delimited by space)\n"
        "    [eg: '0101-XXX-M3-A1_1.fq.gz' prefix is '0101-XXX-M3-A1']",
        required=True,
        nargs="+",
    )
    required.add_argument(
        "-outdir",
        metavar="PATH",
        help="    The path to directory of all output files",
        required=True,
    )
    required.add_argument(
        "-config",
        metavar="STR",
        help="    Config file. Some parameters and paths that need to be set",
        required=True,
    )
    required.add_argument(
        "-bed",
        metavar="STR",
        help="    Input bed file of primers",
        required=True,
    )
    required.add_argument(
        "-baseline-type",
        metavar="STR",
        help="    Select output baseline type. Possible values:{CNV, MSI}",
        choices=["CNV", "MSI"],
        required=True,
    )
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "-baseline-name",
        metavar="STR",
        help="    The name of the CNV/MSI baseline \n"
        "    [default: XXXBaseline_XXXpanel_AtomSeq.XXXgene]",
        default="XXXBaseline_XXXpanel_AtomSeq.XXXgene",
    )
    optional.add_argument(
        "-threads",
        metavar="INT",
        help="    Number of threads to use [default: 8]",
        type=int,
        default=8,
    )
    optional.add_argument(
        "-auto-remove",
        action="store_true",
        help="    Remove all temporary files (will remove Mid-files)",
    )
    argv = vars(parser.parse_args())
    return argv


def GetConfigFile(config, mount_paras, image):
    config_dic = {}
    pydir = os.path.dirname(os.path.realpath(__file__))
    with open(config, "r") as f:
        for line in f:
            if re.match("#", line):
                continue
            if re.search("=", line):
                ls = line.strip().split("=")
                config_dic[ls[0]] = ls[1]
                if ls[0] == "Hg38_Fasta_Path":
                    # tmp_path = "/".join(ls[1].split("/")[:-1])
                    tmp_path = os.path.dirname(ls[1])
                    mount_paras += f" -B {tmp_path}"
                    bwa_dir = Path(f"{tmp_path}/Bwa_index")
                    bismark_dir = Path(f"{tmp_path}/Bisulfite_Genome")
                    if not bwa_dir.exists() or not bismark_dir.exists():
                        if not bwa_dir.exists():
                            bwa_dir.mkdir(parents=True)
                        cmd = (
                            f"python3 {pydir}/bin/make_index.py "
                            f"-image {image} -fasta {ls[1]}"
                        )
                        print(
                            f"   {cmd}\n\033[33mNotice\033[0m: Detected that hg38 "
                            "alignment index doesn't exist. Building for you, "
                            "taking approximately 100 minutes "
                        )
                        os.system(cmd)
                else:
                    mount_paras += f" -B {ls[1]}"
    return config_dic, mount_paras


def FastqDispose(fqdir, sample_name):
    sample_info = {}
    for fq in glob.glob(f"{fqdir}/*_*1.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        suffix = "_" + "_".join(fq.split("/")[-1].split("_")[1:])
        if sample not in sample_info:
            sample_info[sample] = {}
            sample_info[sample]["directory"] = fqdir
            sample_info[sample]["prefix"] = sample
            sample_info[sample]["R1_suffix"] = suffix
        else:
            print(
                f"The {sample} sample name in the rawdata directory is duplicated, please check"
            )
            continue
    for fq in glob.glob(f"{fqdir}/*_*2.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        suffix = "_" + "_".join(fq.split("/")[-1].split("_")[1:])
        if sample not in sample_info:
            print(
                f"In the rawdata directory, {sample} only has fq2 but not fq1. Please check"
            )
            exit()
        else:
            sample_info[sample]["R2_suffix"] = suffix
    non_saved = []
    if sample_name:
        for name in sample_name:
            if name not in list(sample_info.keys()):
                non_saved.append(name)
                continue
        if len(non_saved) > 0:
            print(
                f"Traceback: fastq directory don`t have {'、'.join(non_saved)}, please check"
            )
            exit()
        for sample in list(sample_info.keys()):
            if sample not in sample_name:
                del sample_info[sample]
                continue
    # 检出字典是否为空
    if not sample_info:
        print(
            f"{fqdir} fastq format has wrong or {'、'.join(sample_name)} has wrong, please check"
        )
        exit()
    return sample_info


def defaultYaml(
    sample_dic,
    config_dic,
    outd,
    image,
    baselinetype,
    baselinename,
    autoremove,
    primerbed,
):
    fasta_dir = os.path.dirname(f"{config_dic['Hg38_Fasta_Path']}")
    default = {
        "samples": sample_dic,
        "python3": "python3",
        "scripts": "/snakemake/bin/dev",
        "docs": "/snakemake/resources",
        # "scripts": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/bin/dev",
        # "docs": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/resources",
        # "bed": primerbed,
        "image": image,
        "bed": primerbed,
        "baseline-type": baselinetype,
        "baseline-name": baselinename,
        "auto-remove": autoremove,
        "bedtools": "bedtools",
        "fastp": {
            "path": "fastp",
            "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "threads": 8,
        },
        "bwa": {
            "path": "bwa",
            "fasta": config_dic["Hg38_Fasta_Path"],
            "bwa_index": f"{fasta_dir}/Bwa_index/hg38",
            # "bwa_index": config_dic["Bwa_Index_Path"],
            "read_group": r' "@RG\tID:flowcell1.lane1\tLB:library1\tSM:sample\tPL:ILLUMINA" ',
            "threads": 8,
        },
        "samtools": {"path": "samtools", "threads": 8},
        "gatk": {
            "path": "gatk",
            "known_snps": f"{config_dic['Variant_library']}/GATK_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "known_indels": f"{config_dic['Variant_library']}/GATK_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        },
        "gencore": {"path": "gencore", "support_reads": 1},
        "cnvkit": {
            "path": "cnvkit.py",
            "Rscript_path": "/snakemake/bin/dev/Rscript",
            # "Rscript_path": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/bin/dev/Rscript",
        },
        "msisensor_pro": {"path": "msisensor-pro"},
    }
    with open(f"{outd}/Baseline_config.yaml", "w") as fi:
        yaml.dump(
            default, fi, default_flow_style=False, allow_unicode=True, sort_keys=False
        )
    return default


def turn_primer_to_region(primerbed, outd, mount_paras, image):
    with open(primerbed) as f, open(f"{outd}/zz.region.xls", "w") as fi:
        for line in f:
            ls = line.strip().split("\t")
            if ls[5] == "+":
                fi.write("{}\t{}\t{}\n".format(ls[0], ls[2], str(int(ls[2]) + 80)))
            elif ls[5] == "-":
                fi.write("{}\t{}\t{}\n".format(ls[0], str(int(ls[1]) - 80), ls[1]))
            else:
                print("Wrong: BED format has wrong, please check")
                exit()
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"bedtools sort -i {outd}/zz.region.xls > {outd}/zz.region.bed1"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"bedtools merge -i {outd}/zz.region.bed1 > {outd}/zz.region.bed"
    )
    # print(cmd)
    os.system(cmd)


def pipelineanalysis(
    config_dic,
    outd,
    image,
    baselinetype,
    baselinename,
    autoremove,
    threads,
    mount_paras,
    pydir,
    primerbed,
    fq_lis,
):
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"snakemake -s {pydir}/bin/snakefile_preprocess "
        f"--configfile {outd}/Baseline_config.yaml --cores {threads} --rerun-incomplete "
        f"-d {outd}/Mid-files -q progress --use-conda --printshellcmds"
    )
    # print(cmd)
    os.system(cmd)
    bam_str = ""
    for fq in fq_lis:
        bam_str += f"{outd}/Mid-files/{fq}/{fq}_deduplicated.bam "
    turn_primer_to_region(primerbed, outd, mount_paras, image)  # method
    if baselinetype == "CNV":
        cmd = (
            f"singularity exec -B {mount_paras} {image} "
            f"cnvkit.py target {outd}/zz.region.bed --annotate {pydir}/docs/refFlat.txt "
            f"-o {outd}/zz.tmp.bait 2> {outd}/zz.log"
        )
        # print(cmd)
        os.system(cmd)
        cmd = (
            f"python3 {pydir}/bin/mark_bed_cnvfusion.py --input {outd}/zz.tmp.bait --cnv-gene "
            f"{pydir}/docs/general_CNVgene.txt --output {outd}/zz.tmp.bait.bed"
        )
        # print(cmd)
        os.system(cmd)
        cmd = (
            f"singularity exec -B {mount_paras} {image} "
            f"cnvkit.py batch -n {bam_str} -m amplicon -f {config_dic['Hg38_Fasta_Path']} "
            f"-t {outd}/zz.tmp.bait.bed --target-avg-size 50 "
            f"-d {outd}/Mid-files/CNVcalling --output-reference {outd}/{baselinename} "
            f"2> {outd}/zz.log"
        )
        # print(cmd)
        os.system(cmd)
    else:
        with open(f"{outd}/zz.configure.xls", "w") as fi:
            for fq in fq_lis:
                fi.write(
                    "{}\t{}\n".format(
                        fq, f"{outd}/Mid-files/{fq}/{fq}_deduplicated.bam"
                    )
                )
        cmd = (
            f"""awk -F "\\t" -v OFS="\\t" 'NR>1{{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' """
            f"{pydir}/docs/GRCh38.microsatellites.tsv >{outd}/zz.msi.tsv.nohead1"
        )
        # print(cmd)
        os.system(cmd)
        cmd = (
            f"singularity exec -B {mount_paras} {image} "
            f"bedtools intersect -a {outd}/zz.msi.tsv.nohead1 -b "
            f"{outd}/zz.region.bed -wa |cut -f 1,3- >{outd}/zz.msi.tsv.nohead"
        )
        # print(cmd)
        os.system(cmd)
        cmd = f"head -n 1 {pydir}/docs/GRCh38.microsatellites.tsv|cat - {outd}/zz.msi.tsv.nohead >{outd}/zz.msi.tsv"
        # print(cmd)
        os.system(cmd)
        cmd = (
            f"singularity exec -B {mount_paras} {image} "
            f"msisensor-pro baseline -d {outd}/zz.msi.tsv "
            f"-i {outd}/zz.configure.xls -o {outd}/Mid-files/MSIcalling -c 15 1> {outd}/zz.log"
        )
        # print(cmd)
        os.system(cmd)
        cmd = (
            f"cp  {outd}/Mid-files/MSIcalling/zz.msi.tsv_baseline {outd}/{baselinename}"
        )
        # print(cmd)
        os.system(cmd)
    os.system(f"rm -f {outd}/zz.* ")
    if autoremove:
        midfile = Path(f"{outd}/Mid-files")
        if midfile.exists():
            os.system(f"rm -rf {outd}/Mid-files")


def func(
    image,
    fqdir,
    fq_lis,
    outd,
    configfile,
    primerbed,
    baselinetype,
    baselinename,
    threads,
    autoremove,
):
    pydir = os.path.dirname(os.path.realpath(__file__))
    if not image:
        image = glob.glob(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "*.sif")
        )[0]
    else:
        image = os.path.abspath(image)

    # check abspath
    softlink, hardlink = [], []
    for fq in glob.glob(f"{fqdir}/*_*1.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        if os.path.islink(fq):
            softlink_path = os.path.dirname(fq)
            original_path = os.path.dirname(os.readlink(fq))
            if sample in fq_lis:
                if softlink_path not in softlink:
                    softlink.append(softlink_path)
                if original_path not in hardlink:
                    hardlink.append(original_path)
        else:
            original_path = os.path.dirname(fq)
            if sample in fq_lis:
                if original_path not in hardlink:
                    hardlink.append(original_path)
    mount_path = hardlink.copy()
    outd = os.path.abspath(outd)
    configfile = os.path.abspath(configfile)
    primerbed = os.path.abspath(primerbed)
    outputdir = Path(outd)
    if not outputdir.exists():
        outputdir.mkdir(parents=True)
    mount_path += [outd, configfile, primerbed]
    if softlink:
        mount_path += softlink
    mount_paras = " -B ".join(mount_path)

    config_dic, mount_paras = GetConfigFile(configfile, mount_paras, image)  # method
    sample_dic = FastqDispose(fqdir, fq_lis)  # method
    defaultDic = defaultYaml(
        sample_dic,
        config_dic,
        outd,
        image,
        baselinetype,
        baselinename,
        autoremove,
        primerbed,
    )  # method
    pipelineanalysis(
        config_dic,
        outd,
        image,
        baselinetype,
        baselinename,
        autoremove,
        threads,
        mount_paras,
        pydir,
        primerbed,
        fq_lis,
    )  # method
    time_end = datetime.datetime.now()
    print(
        f"Total time spent: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    )


def main():
    argv = argparse_line()
    func(
        argv["image"],
        argv["fq_dir"],
        argv["fq_prefix"],
        argv["outdir"],
        argv["config"],
        argv["bed"],
        argv["baseline_type"],
        argv["baseline_name"],
        argv["threads"],
        argv["auto_remove"],
    )


if __name__ == "__main__":
    main()
