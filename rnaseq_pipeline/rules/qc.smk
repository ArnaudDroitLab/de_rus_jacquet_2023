
rule fastqc:
    input:
        "{path}/results/{dir}/{sample_pair}.fastq"
    output:
        html = "{path}/results/qc/{dir}/fastqc/{sample_pair}_fastqc.html",
        zip = "{path}/results/qc/{dir}/fastqc/{sample_pair}_fastqc.zip"
    params:
        tmp = expand("{{path}}/{tmp}", tmp = config['qc']['tmp']),
        dir = "{path}/results/qc/{dir}/fastqc"
    log:
        out = "{path}/logs/qc/{dir}/fastqc/{sample_pair}.log",
        err = "{path}/logs/qc/{dir}/fastqc/{sample_pair}.err"
    resources:
        time=lambda wildcards, attempt: attempt * config['qc']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['qc']['mem']
    threads:
        config['qc']['threads']
    container:
        config['qc']['container']
    shell:
        """
        if [ ! -d {params.tmp} ]; 
        then 
            mkdir -p {params.tmp};
        fi; 
        fastqc  -t {threads} \
                {input} \
                -d {params.tmp} \
                -o {params.dir} \
                > {log.out} \
                2> {log.err} \
         """


rule multiqc:
    input:
        zip = expand("{{path}}/results/qc/{{dir}}/fastqc/{sample_pair}_fastqc.zip", 
                        sample_pair = config['rq'].keys())
    output:
        report("{path}/results/qc/{dir}/multiqc/multiqc_{dir}.html", caption="../report/qc.rst",
                    category="QC Step", subcategory="{dir}")
    params:
        outdir = "{path}/results/qc/{dir}"
    log:
        out = "{path}/logs/qc/{dir}/multiqc.log",
        err = "{path}/logs/qc/{dir}/multiqc.err"
    resources:
        time=lambda wildcards, attempt: attempt * config['qc']['time'],
        mem=lambda wildcards, attempt: attempt * 2**10 * config['qc']['mem']
    container:
        config['qc']['container']
    shell:
        """
        multiqc -p -f \
                -n {params.outdir}/multiqc/multiqc_{wildcards.dir} \
                {params.outdir} \
                > {log.out} \
                2> {log.err} \
        """
