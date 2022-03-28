
def get_input(wildcards):
    sample = wildcards.sample
    path = wildcards.path
    pairing = units.loc[sample, 'fq']
    r1 = f"{path}/results/merge/{sample}_R1.fastq"
    r2 = f"{path}/results/merge/{sample}_R2.fastq"
    if pairing == 'PE':
        return r1, r2
    elif pairing == 'SE':
        return r1 

## Not useful right now because of touch output. 
## Could be used to validate output with shell loop
## i.e. : for file in $output; do if [ ! -e $file ]; ...
def get_output(wildcards):
    sample = wildcards.sample
    path = wildcards.path
    pairing = units.loc[sample, 'fq']
    r1 = f"{path}/results/trimming/{sample}_R1.fastq"
    r2 = "{path}/results/trimming/{sample}_R2.fastq"
    r1_unp = "{path}/results/trimming/unpaired/{sample}_R1.fastq"
    r2_unp = "{path}/results/trimming/unpaired/{sample}_R2.fastq"
    if pairing == 'PE':
        return r1, r2, r1_unp, r2_unp

    elif pairing == 'SE':
        return r1, r1_unp


def get_params(wildcards):
    sample = wildcards.sample
    path = wildcards.path
    pairing = units.loc[sample, 'fq']
    r1 = f"{path}/results/merge/{sample}_R1.fastq"
    r2 = f"{path}/results/merge/{sample}_R2.fastq"
    or1 = f"{path}/results/trimming/{sample}_R1.fastq"
    or2 = f"{path}/results/trimming/{sample}_R2.fastq"
    or1_unp = f"{path}/results/trimming/unpaired/{sample}_R1.fastq"
    or2_unp = f"{path}/results/trimming/unpaired/{sample}_R2.fastq"
    if pairing == 'PE':
        return f"-i {r1} -I {r2} -o {or1} -O {or2} --unpaired1 {or1_unp} --unpaired2 {or2_unp}"

    elif pairing == 'SE':
        return f"-i {r1} -o {or1} --unpaired1 {or1_unp}"

rule trimming_fastp:
    input:
        get_input
    output:
        touch("{path}/results/trimming/{sample}.trimmed")
    params:
        io = get_params,
    log:
        json = "{path}/results/trimming/reports/{sample}.json", 
        html = "{path}/results/trimming/reports/{sample}.html",
        out  = "{path}/logs/trimming/{sample}.out",
        err  = "{path}/logs/trimming/{sample}.err",
    threads: 
        config['fastp']['threads']
    resources:
        time=lambda wildcards, attempt: attempt * config['fastp']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['fastp']['mem']
    container: 
        config['fastp']['container']
    shell:
        """
        fastp   {params.io} \
                -h {log.html} \
                -j {log.json} \
                -w {threads} \
                > {log.out} \
                2> {log.err}
        """


