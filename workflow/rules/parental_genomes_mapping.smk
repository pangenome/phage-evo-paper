# [[file:../../main.org::*Map reads to minia assembly][Map reads to minia assembly:1]]
rule prefix_fastq:
    input:
        samples=expand(join_path(config['data']['reads'], '{sample}.merged.fastq'), sample=SAMPLES),
    params:
        samples_prefixed=join_path(config['data']['reads'], 'P1-10.merged.prefixed.fastq')
    output:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.fastq.gz')
    threads:
        get_cores_perc(1)
    shell:
        """
        echo {input.samples} \
            | tr ' ' '\\n' \
            | while read sample; do
                prefix=$( basename $sample | cut -d'.' -f1)
                sed -r '/^@.+runid/ s/^@/@'$prefix'#1#/' $sample >> {params.samples_prefixed}
            done
        cat {params.samples_prefixed} | pigz -p {threads} > {output.samples_prefixed_gzipped}
        """
# Map reads to minia assembly:1 ends here

# [[file:../../main.org::*Filter out bacterial reads][Filter out bacterial reads:1]]
rule minimap2_map_reads_to_refs:
    input:
        target = config['data']['genomes']['merged'],
        samples_prefixed_gzipped = join_path(config['data']['reads'], 'P1-10.merged.prefixed.fastq.gz'),
    output:
        all_reads_sam = join_path(config['data']['reads'], 'P1-10.merged.prefixed_X_genomes.sam'),
    threads:
        get_cores_perc(1)
    conda:
        '../envs/minimap2_env.yaml'
    shell:
        "minimap2 -ax map-pb -t {threads} {input.target} {input.samples_prefixed_gzipped} > {output.all_reads_sam}"
# Filter out bacterial reads:1 ends here
