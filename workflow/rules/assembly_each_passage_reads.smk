# [[file:../../main.org::*Prefix reads][Prefix reads:1]]
rule prefix_fastq:
    input:
        sample = join_path(config['data']['reads'], '{sample}.merged.fastq'),
    output:
        sample_prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz')
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "prefix=$( basename {input.sample} | cut -d'.' -f1) && "
        "sed -r '/^@.+runid/ s/^@/@'$prefix'#1#/' {input.sample} | bgzip > {output.sample_prefixed}"
# Prefix reads:1 ends here

# [[file:../../main.org::*nanoplot][nanoplot:1]]
rule quality_check_plot_before_filtering:
    input:
        sample = join_path(config['data']['reads'], '{state}', '{sample}.{state}.fastq.gz')
    output:
        plot_dir = directory("results/single/nanoplot/{state}/{sample}")
    threads:
        get_cores_perc(1)
    conda:
        "../envs/nanoplot_env.yaml"
    shell:
        "NanoPlot -t 4 --plots dot -o {output.plot_dir} --fastq {input.sample}"
# nanoplot:1 ends here

# [[file:../../main.org::*FILTER READS][FILTER READS:1]]
rule filter_reads:
    input:
        prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz')
    output:
        filtered = join_path(config['data']['reads'], 'filtered', '{sample}.filtered.fastq.gz')
    params:
        **config['params']['filtlong']
    conda:
        "../envs/filtlong_env.yaml"
    threads:
        10
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input.prefixed} | pigz > {output.filtered}"
# FILTER READS:1 ends here

# [[file:../../main.org::*Minia assembly][Minia assembly:1]]
rule minia:
    input:
        filtered = join_path(config['data']['reads'], 'filtered', '{sample}.filtered.fastq.gz')
    output:
        minia_assembly =  join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.fa")
    threads:
        4
    params:
        **config['params']['minia'],
        prefix_fasta=minia_prefix
    conda:
        '../envs/minia_env.yaml'
    shell:
        "minia -nb-cores {params.cores} -kmer-size {params.kmer} -abundance-min {params.abundance} -out {params.prefix_fasta} -in {input.filtered} && "
        " find $( dirname {output.minia_assembly} ) -type f ! -name '*'$(basename {output.minia_assembly}) -exec rm {{}} \;"
# Minia assembly:1 ends here

# [[file:../../main.org::*fasta to gfa][fasta to gfa:1]]
rule minia_fasta_to_gfa:
    input:
        minia_assembly =  join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.fa"),
        script=join_path(snakefile_path, 'scripts', 'convertToGFA.py'),
    output:
        minia_assembly_gfa = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.gfa")
    params:
        **config['params']['minia'],
    conda:
        '../envs/minia_env.yaml'
    threads:
        10
    shell:
        "python {input.script} {input.minia_assembly} {output.minia_assembly_gfa} {params.kmer}"
# fasta to gfa:1 ends here

# [[file:../../main.org::*Graphaligner MINIA][Graphaligner MINIA:1]]
rule polishing_graphaligner_minia:
    input:
        samples_prefixed_gzipped = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz'),
        minia_assembly_gfa = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.gfa")
    output:
        minia_gaf = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.polished.gaf"),
        minia_assembly_gfa_polished = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.polished.fa"),
    threads:
        4
    params:
        dbtype = "vg",
        seed_minimizer = 15
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        "GraphAligner -g {input.minia_assembly_gfa} -f {input.samples_prefixed_gzipped} -x {params.dbtype} --threads 10 --seeds-minimizer-length {params.seed_minimizer} --seeds-minimizer-windowsize {params.seed_minimizer} -a {output.minia_gaf} --corrected-out {output.minia_assembly_gfa_polished}"
# Graphaligner MINIA:1 ends here

# [[file:../../main.org::*Filter by length][Filter by length:1]]
rule filter_by_length:
    input:
        minia_assembly_gfa_polished = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.polished.fa"),
        script = join_path(snakefile_path, 'scripts', 'filter_by_length.py')
    output:
        minia_assembly_polished_filtered = join_path('results', 'single', 'minia', '{sample}', '{sample}.' + minia_prefix + ".contigs.polished" + filter_contigs_prefix + ".fa.gz")
    params:
        **config['params']['minia']
    conda:
        '../envs/bio_env.yaml'
    threads:
        10
    shell:
        "python3 {input.script} {input.minia_assembly_gfa_polished} {params.min_contig_lenght}  {params.max_contig_lenght} | bgzip > {output.minia_assembly_polished_filtered}"
# Filter by length:1 ends here
