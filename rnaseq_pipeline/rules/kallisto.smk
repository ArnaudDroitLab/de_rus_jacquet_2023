
rule kallisto_index:
    input:
        fasta = config["resources"]["ref"]["transcriptome"]
    output:
        index = "{path}/results/kallisto/transcripts.idx"
    params:
        extra = "--make-unique"
    log:
        out = "{path}/logs/kallisto/index.log",
        err = "{path}/logs/kallisto/index.err"
    threads: 1
    resources:
        time=lambda wildcards, attempt: attempt * config['kallisto']['index']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['kallisto']['index']['mem']
    container: 
        config['kallisto']['container']
    shell:
        """
        kallisto    index \
                    {params.extra} \
                    --index={output.index} \
                    {input.fasta} \
                    >  {log.out} \
                    2> {log.err}
        """

def get_kallisto_input(wildcards):
    sample = wildcards.sample
    path = wildcards.path
    pairing = units.loc[sample, 'fq']
    r1 = f"{path}/results/trimming/{sample}_R1.fastq"
    r2 = f"{path}/results/trimming/{sample}_R2.fastq"
    if pairing == 'PE':
        return f"{r1} {r2}"
    elif pairing == 'SE':
        return r1 

def kallisto_params(wildcards):
    extra = config["kallisto"]["params"]
    sample = wildcards.sample
    pairing = units.loc[sample, 'fq']
    if pairing == 'SE':
        extra += " --single -l 200 -s 30"
    # else:
    #     extra += " --fusion"
    return extra

rule kallisto_quant:
    input:
        idx = "{path}/results/kallisto/transcripts.idx",
        flag = "{path}/results/trimming/{sample}.trimmed"
    output:
        expand("{{path}}/results/kallisto/{{sample}}/{file}", file = ["abundance.tsv", "abundance.h5", "run_info.json"])
    log:
        out = "{path}/logs/kallisto/quant/{sample}.log",
        err = "{path}/logs/kallisto/quant/{sample}.err"
    params:
        fq = get_kallisto_input,
        output = "{path}/results/kallisto/{sample}",
        extra = kallisto_params
    container:
        config['kallisto']['container']
    threads: 
        config['kallisto']['quant']['threads']
    resources:
        time=lambda wildcards, attempt: attempt * config['kallisto']['quant']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['kallisto']['quant']['mem']
    shell:
        """
        kallisto    quant \
                    -i {input.idx} \
                    -o {params.output} \
                    -t {threads} \
                    {params.extra} \
                    {params.fq} \
                    >  {log.out} \
                    2> {log.err}
        """


