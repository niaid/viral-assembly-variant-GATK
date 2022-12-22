shell.prefix("set -eo pipefail; ")

wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index)

def myref(param):
    if param == "denovo":
        return "data/{sample}/ragtag/ragtag.scaffold.clean.fasta"
    elif param == "mapping":
        return config["params"]["ref"]

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")

def get_filter_info(wildcards):  
    return [config["filtering"]["hard"][wildcards.vartype]["info"]]

def get_filt_samples(wildcards):
    return expand("variants/filtered/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards)

def get_filt_merge(wildcards):
    return " ".join("INPUT={}".format(f) for f in expand("variants/filtered/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards))

    