
def get_input(wildcards):
    sample = wildcards.sample_pair
    return config['rq'][sample]

def get_commandline(wildcards):
    sample = wildcards.sample_pair
    path = wildcards.path
    
    nb_samples = len(config['rq'][sample])
    samples = ' '.join(config['rq'][sample])
    
    cmd  = f"cat {samples} > {path}/results/merge/{sample}.fastq"
    if nb_samples == 1:
        cmd  = f"ln -s {samples} {path}/results/merge/{sample}.fastq"

    return cmd


rule merge_lanes:
    input:
        get_input
    output:
        "{path}/results/merge/{sample_pair}.fastq"
    params:
        cmd = get_commandline
    resources:
        mem=lambda wildcards, attempt: attempt * 2**10 * 1,
        time=lambda wildcards, attempt: attempt * 30
    shell:
        """
        {params.cmd}
        """
