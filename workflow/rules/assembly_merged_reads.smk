# [[file:../../main.org::*Merge reads][Merge reads:1]]
rule prefix_fastq:
    input:
        samples=expand(join_path(config['data']['reads'], '{sample}.merged.fastq'), sample=SAMPLES),
    params:
        samples_prefixed = join_path(config['data']['reads'], 'P1-10.merged.prefixed.before_qc.fastq'),
    output:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.before_qc.fastq.gz'),
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
        pigz -p {threads} {params.samples_prefixed}
        """
# Merge reads:1 ends here

# [[file:../../main.org::*nanoplot][nanoplot:1]]
rule nanoplot:
    input:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.{state}_qc.fastq.gz'),
    output:
        directory("results/nanoplot/{state}_filter")
    threads:
        get_cores_perc(0.5)
    conda:
        "../envs/nanoplot_env.yaml"
    shell:
        "NanoPlot -t {threads} --plots dot -o {output} --fastq {input}"
# nanoplot:1 ends here

# [[file:../../main.org::*FILTER READS][FILTER READS:1]]
rule filter_reads:
    input:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.before_qc.fastq.gz'),
    output:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.after_qc.fastq.gz'),
    params:
        **config['params']['filtlong']
    conda:
        "../envs/filtlong_env.yaml"
    threads:
        get_cores_perc(0.2)
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input} | pigz -p {threads} > {output}"
# FILTER READS:1 ends here

# [[file:../../main.org::*MINIA3][MINIA3:1]]
rule minia:
    input:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.after_qc.fastq.gz'),
    output:
        minia_assembly=minia_prefix+".contigs.fa"
    threads:
        get_cores_perc(1)
    params:
        **config['params']['minia'],
        prefix_fasta=minia_prefix
    conda:
        '../envs/minia_env.yaml'
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
        '../envs/minia_env.yaml'
    shell:
        "python {input.script} {input.minia_assembly} {output.minia_assembly_gfa} {params.kmer}"
# FASTA_TO_GFA:1 ends here

# [[file:../../main.org::*Graphaligner MINIA][Graphaligner MINIA:1]]
rule polishing_graphaligner_minia:
    input:
        samples_prefixed_gzipped = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz'),
        minia_assembly_gfa = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.gfa')
    output:
        minia_gaf = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished.gaf'),
        minia_assembly_gfa_polished = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished.fa'),
    threads:
        get_cores_perc(1)
    params:
        dbtype = "vg",
        seed_minimizer = 15
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        "GraphAligner -g {input.minia_assembly_gfa} -f {input.samples_prefixed_gzipped} -x {params.dbtype} "
        "--threads {threads} --seeds-minimizer-length {params.seed_minimizer} "
        "--seeds-minimizer-windowsize {params.seed_minimizer} "
        "-a {output.minia_gaf} --corrected-out {output.minia_assembly_gfa_polished}"
# Graphaligner MINIA:1 ends here

# [[file:../../main.org::*Filter by length][Filter by length:1]]
rule filter_by_length_and_index:
    input:
        minia_assembly_gfa_polished = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished.fa'),
        script = join_path(snakefile_path, 'scripts', 'filter_by_length.py')
    output:
        minia_assembly_polished_filtered = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished' + filter_contigs_prefix + ".fa.gz"),
        fai = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished' + filter_contigs_prefix + ".fa.gz.fai"),
        gzi = join_path('results', results_dir, 'minia', '{sample}', '{sample}.contigs.polished' + filter_contigs_prefix + ".fa.gz.gzi"),
    params:
        **config['params']['minia']
    conda:
        '../envs/bio_env.yaml'
    threads:
        1
    shell:
        "python3 {input.script} {input.minia_assembly_gfa_polished} {params.min_contig_lenght}  {params.max_contig_lenght} | bgzip > {output.minia_assembly_polished_filtered} && "
        "samtools faidx {output.minia_assembly_polished_filtered}"
# Filter by length:1 ends here

# [[file:../../main.org::*Create index][Create index:1]]
rule create_index_fasta:
    input:
        minia_assembly_polished_filtered = filter_contigs_prefix + '.contigs.polished.fa',
    output:
        minia_assembly_polished_filtered_crompressed = filter_contigs_prefix + '.contigs.polished.fa.gz',
        fai = filter_contigs_prefix + '.contigs.polished.fa.gz.fai',
        gzi = filter_contigs_prefix + '.contigs.polished.fa.gz.gzi',
    threads:
        get_cores_perc(0.5)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "cat {input.minia_assembly_polished_filtered} | bgzip -@ {threads} > {output.minia_assembly_polished_filtered_crompressed} && "
        "samtools faidx {output.minia_assembly_polished_filtered_crompressed}"
# Create index:1 ends here

# [[file:../../main.org::*Get sample and add parental phages genomes][Get sample and add parental phages genomes:1]]
rule add_parental_genomes_and_get_sample:
    input:
        minia_assembly_polished_filtered_crompressed = filter_contigs_prefix + '.contigs.polished.fa.gz',
        parental_genomes = config['data']['genomes']['ecoli_and_phages']
    params:
        prefix = filter_contigs_prefix + '.contigs.polished.sample1K.fa',
    output:
        minia_assembly_polished_filtered_crompressed_sampled = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz',
        fai = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz.fai',
        gzi = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz.gzi',
    threads:
        get_cores_perc(0.8)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "cat {input.parental_genomes} > {params.prefix} && "
        "samtools faidx {input.minia_assembly_polished_filtered_crompressed} "
        "$( seq 1 10 | while read i; do zgrep  -P '^>P'$i'#' {input.minia_assembly_polished_filtered_crompressed} | shuf -n 100 ; done | sed 's/>//' ) "
        ">> {params.prefix} && "
        " bgzip -@ {threads}  {params.prefix} && "
        " samtools faidx {output.minia_assembly_polished_filtered_crompressed_sampled}"
# Get sample and add parental phages genomes:1 ends here

# [[file:../../main.org::*PGGB minia_polished][PGGB minia_polished:1]]
rule pggb_minia:
    input:
        minia_assembly_polished_filtered_crompressed_sampled = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz',
        fai = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz.fai',
        gzi = filter_contigs_prefix + '.contigs.polished.sample1K.fa.gz.gzi',
    output:
        directory( "results/pggb/minia.assembly" + pggb_prefix + ".ecoli.and.phages" ),
    params:
        **config['params']['pggb']
    conda:
        '../envs/pggb_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        "n_mappings=$( zgrep -c '>' {input.minia_assembly_polished_filtered_crompressed_sampled} ) && "
        " pggb -m -p {params.map_pct_id} -n $n_mappings -s {params.segment_length} -l {params.block_length} -t {threads} -o {output} -i {input.minia_assembly_polished_filtered_crompressed_sampled}"
# PGGB minia_polished:1 ends here

# [[file:../../main.org::*Get distance matrix][Get distance matrix:1]]
rule odgi_get_distance_matrix:
    input:
        odgi_graph = glob.glob(join_path("results/pggb/minia.assembly" + pggb_prefix + ".ecoli.and.phages",  '*.smooth.final.og'))[0]
    output:
# Get distance matrix:1 ends here
