# [[file:../../main.org::*config][config:1]]
results_dir = os.path.basename(workflow.snakefile).replace('.smk', '')
# config:1 ends here

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

# [[file:../../main.org::*FASTQ lenght filtering][FASTQ lenght filtering:1]]
rule filter_reads:
    input:
        prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz')
    output:
        filtered = join_path(config['data']['reads'], 'filtered', '{sample}.filtered.fastq.gz')
    params:
        **config['params']['seqkit']
    conda:
        "../envs/seqkit_env.yaml"
    threads:
        4
    shell:
        "seqkit seq {input.prefixed} -j {threads} -m {params.min} -M {params.max} | bgzip > {output.filtered}"
# FASTQ lenght filtering:1 ends here

# [[file:../../main.org::*Minia assembly][Minia assembly:1]]
rule minia:
    input:
        prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz'),
        script_abundance = join_path(snakefile_path, 'scripts', 'get_abundance.sh'),
    output:
        minia_assembly =  join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.fa'),
    threads:
        4
    params:
        **config['params']['minia'],
    conda:
        '../envs/minia_env.yaml'
    shell:
        "RELATIVE_ABUNDACE=$( {input.script_abundance} {params.P1_abundance} {params.P1_bp} {input.prefixed} ) && "
        "minia -nb-cores {params.cores} -kmer-size {params.kmer} -abundance-min $RELATIVE_ABUNDACE "
        "-out $(echo {output.minia_assembly} | sed 's/.contigs.fa//') -in {input.prefixed} && "
        "find $( dirname {output.minia_assembly} ) -type f ! -name '*'$(basename {output.minia_assembly}) -exec rm {{}} \;"
# Minia assembly:1 ends here

# [[file:../../main.org::*fasta to gfa][fasta to gfa:1]]
rule minia_fasta_to_gfa:
    input:
        minia_assembly =  join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.fa'),
        script = join_path(snakefile_path, 'scripts', 'convertToGFA.py'),
    output:
        minia_assembly_gfa = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.gfa')
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
        filtered = join_path(config['data']['reads'], 'filtered', '{sample}.filtered.fastq.gz'),
        minia_assembly_gfa = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.gfa'),
    output:
        minia_gaf = join_path('results', results_dir, 'minia', '{sample}', '{sample}.reads.polished.gaf'),
        polished_reads = join_path('results', results_dir, 'minia', '{sample}', '{sample}.reads.polished.fa'),
    threads:
        get_cores_perc(0.3)
    params:
        dbtype = "vg",
        seed_minimizer = 15
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        "GraphAligner -g {input.minia_assembly_gfa} -f {input.filtered} -x {params.dbtype} "
        "--threads {threads} --seeds-minimizer-length {params.seed_minimizer} "
        "--seeds-minimizer-windowsize {params.seed_minimizer} -a {output.minia_gaf} "
        "--corrected-out {output.polished_reads}"
# Graphaligner MINIA:1 ends here

# [[file:../../main.org::*Sample 1000][Sample 1000:1]]
rule sample_genomes:
    input:
        polished_reads = join_path('results', results_dir, 'minia', '{sample}', '{sample}.reads.polished.fa'),
    output:
        polished_reads = join_path('results', results_dir, 'minia', '{sample}', '{sample}.reads.polished.sample.fa.gz' ),
    params:
        sample_size = 100
    threads:
        4
    shell:
        "bgzip -@ {threads} {input.polished_reads} && "
        "samtools faidx {input.polished_reads}.gz $(zgrep '>' {input.polished_reads}.gz | sed 's/>//' | shuf -n {params.sample_size}) | "
        "bgzip > {output.polished_reads}"
# Sample 1000:1 ends here

# [[file:../../main.org::*Merge samples][Merge samples:1]]
rule merge_samples_and_parental_genomes:
    input:
        polished_reads = expand(join_path('results', results_dir, 'minia', '{sample}', '{sample}.reads.polished.sample.fa.gz' ), sample=SAMPLES),
        ecoli_and_phages = config['data']['genomes']['ecoli_and_phages'],
    output:
        pggb_input = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz'),
        fai = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz.fai'),
        gzi = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz.gzi'),
    conda:
        '../envs/pggb_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        "cat {input.ecoli_and_phages} <(zcat {input.polished_reads}) | bgzip -@ {threads} > {output.pggb_input} && "
        "samtools faidx {output.pggb_input}"
# Merge samples:1 ends here

# [[file:../../main.org::*Pangenome PGGB][Pangenome PGGB:1]]
rule pggb_pangenome:
    input:
        pggb_input = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz'),
        fai = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz.fai'),
        gzi = join_path('results', results_dir, 'pggb', 'minia.merged.1K.sample.fa.gz.gzi'),
    output:
        pggb_out = directory(join_path('results', results_dir, 'pggb', 'out')),
    params:
        **config['params']['pggb']
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "n_mappings=$( zgrep -c '>' {input.pggb_input} ) && "
        "pggb -m -p {params.map_pct_id} -n $n_mappings -s {params.segment_length} -l {params.block_length} -t {threads} -o {output.pggb_out} -i {input.pggb_input}"
# Pangenome PGGB:1 ends here

# [[file:../../main.org::*Get distance][Get distance:1]]
rule get_distance_metrics:
    input:
        pggb_out = join_path('results', results_dir, 'pggb', 'out'),
    output:
        distance_tsv = join_path('results', results_dir, 'pggb', 'distance_matrix.tsv'),
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "odgi paths -t {threads} -d -i {input.pggb_out}/*.smooth.final.og > {output.distance_tsv}"
# Get distance:1 ends here

# [[file:../../main.org::*Plot phylogeny][Plot phylogeny:1]]
rule plot_phylogeny:
    input:
        distance_tsv = join_path('results', results_dir, 'pggb', 'distance_matrix.tsv'),
        script_phylogeny = join_path(snakefile_path, 'scripts', 'phylogeny.R'),
    output:
        rectangular = join_path('results', results_dir, 'plots', 'ggtree.ecoli.phages.passages.rectangular.pdf'),
        daylight = join_path('results', results_dir, 'plots', 'ggtree.ecoli.phages.passages.daylight.pdf'),
    conda:
        '../envs/R_env.yaml'
    threads:
        1
    shell:
        "Rscript {input.script_phylogeny} {input.distance_tsv} {output.rectangular}"
# Plot phylogeny:1 ends here