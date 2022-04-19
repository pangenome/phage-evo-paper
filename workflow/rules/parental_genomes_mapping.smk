# [[file:../../main.org::*Map reads to minia assembly][Map reads to minia assembly:1]]
rule map_minia_assembly_on_bacterial:
    input:
        minia_assembly_polished_filtered = filter_contigs_prefix + '.polished.prefixed.fa',
        bacterial_genome = config['data']['genomes']['ecoli']
    output:
        minimap2_bam = join_path('results', 'map_reads', 'minia_assembly_X_bacterial_genome.bam')
    conda:
        '../envs/miniasm_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.bacterial_genome} {input.minia_assembly_polished_filtered} \
            | samtools view -b - > {output.minimap2_bam}
        """
# Map reads to minia assembly:1 ends here
