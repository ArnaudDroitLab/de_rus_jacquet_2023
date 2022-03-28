
rule r_rnaseq:
    input:
        expand("{path}/results/kallisto/{sample}", path = config['path'], sample = units.index),
    output:
        "{path}/results/rnaseq/tpm.csv",
        "{path}/results/rnaseq/txi.rds"
    log:
        out = "{path}/logs/rnaseq.log",
        err = "{path}/logs/rnaseq.err"
    params:
        dir = "{path}/results/kallisto",
        prefix = "{path}/results/rnaseq",
        script = config['rnaseq']['script'],
        anno = config['rnaseq']['anno']
    container:
        config['rnaseq']['container']
    resources:
        time=lambda wildcards, attempt: attempt * config['rnaseq']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['rnaseq']['mem']
    shell:
        """
        Rscript {params.script} \
                {params.dir} \
                {params.prefix} \
                {params.anno} \
                > {log.out} \
                2> {log.err}
        """
