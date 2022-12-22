shell.prefix("set -eo pipefail; ")


accession=config["params"]["ref"]
ref_name=config["params"]["ref_name"]


def get_filt_samples(wildcards):
    return expand("data/{sample}/gatk/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards)

def get_filt_merge(wildcards):
    return " ".join("INPUT={}".format(f) for f in expand("data/{sample}/gatk/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards))

rule index_ref:
    input:
        ref=config["params"]["ref"]
    output:
        multiext(f"{accession}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    log:
        "log/reference/bowtie2_index.log"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        bowtie2-build {input.ref} {input.ref}
    """

rule index_ref_gatk:
    input:
        ref=config["params"]["ref"]
    output:
        faidx=f"{accession}.fai", 
        genome_dict=f"ref/{ref_name}.dict"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        samtools faidx {input.ref}
        samtools dict {input.ref} > {output.genome_dict}
    """

rule bowtie2:
    input:
        in1="data/{sample}/ec/{sample}_R1.clean.ec.fastq.gz",
        in2="data/{sample}/ec/{sample}_R2.clean.ec.fastq.gz",
        ref=config["params"]["ref"],
        idx=rules.index_ref.output

    output:
        "data/{sample}/bowtie2/{sample}.sort.bam"
    params:
        sample="{sample}"
    conda:
        "../envs/consiter.yaml"
    shell: """
    bowtie2 --no-unal \
        --rg-id {params.sample} --rg SM:{params.sample} --rg LB:1 --rg PU:1 --rg PL:Illumina \
        -x {input.ref} -1 {input.in1} -2 {input.in2} | samtools sort - > {output}
    """

rule dedup_bam:
    input:
        "data/{sample}/bowtie2/{sample}.sort.bam"
    output:
        "data/{sample}/bowtie2/{sample}.sort.dedup.bam"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        samtools rmdup {input} {output}
        """

rule index_bam:
    input:
        "data/{sample}/bowtie2/{sample}.sort.dedup.bam"
    output:
        "data/{sample}/bowtie2/{sample}.sort.dedup.bam.bai"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        samtools index {input}
        """

rule call_variants:
    input:
        bam="data/{sample}/bowtie2/{sample}.sort.dedup.bam",
        bai="data/{sample}/bowtie2/{sample}.sort.dedup.bam.bai",
        faidx=f"{accession}.fai",
        ref=config["params"]["ref"]
    output:
        vcf="data/{sample}/gatk/{sample}.vcf.gz"
    params:
        ploidy=config["params"]["ploidy"]
    conda:
        "../envs/consiter.yaml"
    shell:"""
        gatk --java-options "-Xmx4g" HaplotypeCaller  \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -ploidy {params.ploidy}
            """

rule select_calls:
    input:
        ref=config["params"]["ref"],
        vcf="data/{sample}/gatk/{sample}.vcf.gz"
    output:
        vcf=temp("data/{sample}/gatk/{sample}.{vartype}.vcf.gz"),
        tbi=temp("data/{sample}/gatk/{sample}.{vartype}.vcf.gz.tbi")
    params:
        extra=get_vartype_arg
    log:
        "log/gatk/selectvariants/{sample}.{vartype}.log"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        gatk SelectVariants -R {input.ref} -V {input.vcf} {params.extra} -O {output.vcf}
    """

rule hard_filter_calls:
    input:
        ref=config["params"]["ref"],
        vcf="data/{sample}/gatk/{sample}.{vartype}.vcf.gz",
        tbi="data/{sample}/gatk/{sample}.{vartype}.vcf.gz.tbi"
    output:
        vcf=temp("data/{sample}/gatk/{sample}.{vartype}.filt.vcf.gz"),
        tbi=temp("data/{sample}/gatk/{sample}.{vartype}.filt.vcf.gz.tbi")
    params:
        filters_info=get_filter_info
    log:
        "log/gatk/variantfiltration/{sample}.{vartype}.log"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        gatk VariantFiltration \
           -R {input.ref} \
           -V {input.vcf} \
           -O {output.vcf} \
           --filter-name "snv-hard-filter" \
           --filter-expression "{params.filters_info}"
        """

rule merge_ind_calls:
    input:
        vcf=get_filt_samples
    output:
        vcf="data/{sample}/{sample}.gatk.filt.vcf.gz",
        tbi=temp("data/{sample}/{sample}.gatk.filt.vcf.gz.tbi")
    params:
        inputs=get_filt_merge
    log:
        "log/picard/merge-filtered.{sample}.log"
    conda:
        "../envs/consiter.yaml"
    shell:"""
        picard MergeVcfs \
            {params.inputs} \
            O={output.vcf}
    """
