# [[file:../../main.org::*ALL][ALL:1]]
rule all:
    input:
        minimap2_bam = join_path('results', 'map_reads', 'minimap2_bacterial_and_phage_genomes.bam')
# ALL:1 ends here

# [[file:../../main.org::*Merge reads][Merge reads:1]]
rule merged_reads:
    input:
        expand(config['reads']+'/{sample}.fastq', sample=SAMPLES)
    output:
        config['data']+'/merged.before.filtering.fastq'
    shell:
        "cat {input} > {output}"
# Merge reads:1 ends here

# [[file:../../main.org::*NANOPLOT][NANOPLOT:1]]
rule nanoplot:
    input:
        config['data']+"/merged.{status}.filtering.fastq"
    output:
        directory("outputs/nanoplot/{status}_filter")
    threads:
        get_cores_perc(0.5)
    conda:
        "envs/nanoplot_env.yaml"
    shell:
        "NanoPlot -t {threads} --plots dot -o {output} --fastq {input}"
# NANOPLOT:1 ends here

# [[file:../../main.org::*FILTER READS][FILTER READS:1]]
rule filter_reads:
    input:
        config['data']+'/merged.before.filtering.fastq'
    output:
        config['data']+'/merged.after.filtering.fastq'
    params:
        **config['params']['filtlong']
    conda:
        "envs/filtlong_env.yaml"
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input} > {output} "
# FILTER READS:1 ends here

# [[file:../../main.org::*MINIA3][MINIA3:1]]
rule minia:
    input:
        config['data']+'/merged.after.filtering.fastq'
    output:
        minia_assembly=minia_prefix+".contigs.fa"
    threads:
        get_cores_perc(0.5)
    params:
        **config['params']['minia'],
        prefix_fasta=minia_prefix
    conda:
        'envs/minia_env.yaml'
    shell:
        "minia -nb-cores {threads} -kmer-size {params.kmer} -abundance-min {params.abundance} -out {params.prefix_fasta} -in {input}"
# MINIA3:1 ends here

# [[file:../../main.org::*FASTA_TO_GFA][FASTA_TO_GFA:1]]
rule minia_fasta_to_gfa:
    input:
        minia_assembly=minia_prefix+".contigs.fa",
        script=join_path(snakefile_path, 'scripts', 'convertToGFA.py'),
    output:
        minia_assembly_gfa=minia_prefix+'.contigs.gfa'
    params:
        **config['params']['minia'],
    conda:
        'envs/minia_env.yaml'
    shell:
        "python {input.script} {input.minia_assembly} {output.minia_assembly_gfa} {params.kmer}"
# FASTA_TO_GFA:1 ends here

# [[file:../../main.org::*Graphaligner MINIA][Graphaligner MINIA:1]]
rule polishing_graphaligner_minia:
    conda:
        'envs/graphaligner_env.yaml'
    input:
        raw_reads=config['data']+'/merged.before.filtering.fastq',
        minia_assembly_gfa=minia_prefix+'.contigs.gfa'
    output:
        minia_gaf=minia_prefix+'.contigs.gaf',
        minia_assembly_gfa_polished=minia_prefix+'.contigs.polished.fa'
    threads:
        get_cores_perc(1)
    params:
        dbtype = "vg",
        seed_minimizer = 15
    shell:
        "GraphAligner -g {input.minia_assembly_gfa} -f {input.raw_reads} -x {params.dbtype} --threads {threads} --seeds-minimizer-length {params.seed_minimizer} --seeds-minimizer-windowsize {params.seed_minimizer} -a {output.minia_gaf} --corrected-out {output.minia_assembly_gfa_polished}"
# Graphaligner MINIA:1 ends here

# [[file:../../main.org::*Filter by length][Filter by length:1]]
rule filter_by_length:
    input:
        minia_assembly_gfa_polished=minia_prefix+'.contigs.polished.fa',
        script = join_path(snakefile_path, 'scripts', 'filter_by_length.py')
    output:
        minia_assembly_polished_filtered = minia_prefix + '.contigs.polished' + filter_contigs_prefix + '.fa'
    params:
        **config['params']['minia']
    conda:
        'envs/bio_env.yaml'
    shell:
        "python3 {input.script} {params.min_contig_lenght}  {params.max_contig_lenght} > {output.minia_assembly_polished_filtered}"
# Filter by length:1 ends here

# [[file:../../main.org::*PGGB minia_polished][PGGB minia_polished:1]]
rule pggb_minia:
    input:
        minia_assembly_polished_filtered = minia_prefix + '.contigs.polished' + filter_contigs_prefix + '.fa'
    output:
        directory("outputs/pggb/minia"+pggb_prefix)
    params:
        **config['params']['pggb']
    conda:
        'envs/pggb_env.yaml'
    threads:
        get_cores_perc(0.5)
    shell:
        "pggb -i {input.minia_assembly_polished_filtered} -p {params.map-pct-id} -n {params.n-mappings} -s {params.segment-length} -t {threads} -o {output} -m"
# PGGB minia_polished:1 ends here
