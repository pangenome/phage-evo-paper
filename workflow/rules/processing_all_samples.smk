# [[file:../../main.org::*Plot lengths][Plot lengths:1]]
rule quality_check_plot:
    input:
        reads = join_path(config['data']['reads'], '{sample}' + config['reads_suffix'])
    output:
        plot_dir = directory(join_path(results_dir, 'nanoplot', '{sample}'))
    threads:
        4
    conda:
        "../envs/nanoplot_env.yaml"
    shell:
        "NanoPlot -t {threads} --plots dot --fastq {input.reads} -o {output.plot_dir}"
# Plot lengths:1 ends here

# [[file:../../main.org::*Assembly][Assembly:1]]
rule minia_assembly:
    input:
        reads = join_path(config['data']['reads'], '{sample}' + config['reads_suffix']),
        script_abundance = join_path(snakefile_path, 'scripts', 'get_abundance.sh'),
        script_fa_to_gfa = join_path(snakefile_path, 'scripts', 'convertToGFA.py'),
    output:
        minia_assembly = join_path(results_dir, 'minia', '{sample}', '{sample}.minia.contigs.fa'),
        minia_assembly_gfa = join_path(results_dir, 'minia', '{sample}', '{sample}.minia.contigs.gfa'),
        log_abundance = join_path('logs', '{sample}.abundance.txt'),
    params:
        **config['params']['minia'],
    threads:
        get_cores_perc(0.1)
    conda:
        '../envs/minia_env.yaml'
    shell:
        "RELATIVE_ABUNDACE=$( {input.script_abundance} {params.P1_abundance} {params.P1_bp} {input.reads} ) && "
        'echo "{wildcards.sample},${{RELATIVE_ABUNDACE}}" > {output.log_abundance} && '
        "minia -nb-cores {threads} -kmer-size {params.kmer} -abundance-min $RELATIVE_ABUNDACE "
        "-out $(echo {output.minia_assembly} | sed 's/.contigs.fa//') -in {input.reads} && "
        "find $( dirname {output.minia_assembly} ) -type f ! -name '*'$(basename {output.minia_assembly}) -exec rm {{}} \; && "
        "python {input.script_fa_to_gfa} {output.minia_assembly} {output.minia_assembly_gfa} {params.kmer}"
# Assembly:1 ends here

# [[file:../../main.org::*Error correction][Error correction:1]]
rule graphaligner_error_correction:
    input:
        reads = join_path(config['data']['reads'], '{sample}' + config['reads_suffix']),
        minia_assembly_gfa = join_path(results_dir, 'minia', '{sample}', '{sample}.minia.contigs.gfa'),
    output:
        putative_phage_genomes = join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.fastq'),
        putative_phage_genomes_polished = join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.polished' + '.prefixed' + '.fa.gz'),
    params:
        gam = join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.polished' + '.gam'),
        fasta = join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.polished' + '.fa'),
        **config['params']['seqkit'],
        **config['params']['graphaligner'],
    threads:
        get_cores_perc(0.3)
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        "seqkit seq {input.reads} -j {threads} -m {params.min} -M {params.max} > {output.putative_phage_genomes} && "
        "GraphAligner -g {input.minia_assembly_gfa} -f {output.putative_phage_genomes} -x {params.dbtype} "
        "--threads {threads} --seeds-minimizer-length {params.seed_minimizer} "
        "--seeds-minimizer-windowsize {params.seed_minimizer} -a {params.gam} "
        "--corrected-out {params.fasta} && "
        "sed -r '/>/ s|>|>{wildcards.sample}#1#|;s|\s.+||' {params.fasta} | bgzip > {output.putative_phage_genomes_polished} && "
        "rm {params.gam} {params.fasta}"
# Error correction:1 ends here

# [[file:../../main.org::*Sample and merge][Sample and merge:1]]
rule merge_and_sample:
    input:
        putative_phage_genomes_polished = expand(join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.polished' + '.prefixed' + '.fa.gz'), sample=SAMPLES),
    output:
        pggb_input = join_path(results_dir, 'pggb', 'genomes.sample_' + str(config['sample_size']) + '_from_each_passage.fa.gz')
    params:
        fasta = join_path(results_dir, 'pggb', 'genomes.sample_' + str(config['sample_size']) + '_from_each_passage.fa'),
        sample_size = config['sample_size']
    threads:
        get_cores_perc(0.5)
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        '> {params.fasta} && '
        'for f in {input.putative_phage_genomes_polished}; do '
        "samtools faidx $f $( zgrep -Po '(?<=^\>).+' $f | shuf -n {params.sample_size} ) >> {params.fasta}; done && "
        'bgzip -@ {threads} {params.fasta} && samtools faidx {output.pggb_input} '
# Sample and merge:1 ends here

# [[file:../../main.org::*Fastani][Fastani:1]]
rule fastaANI_distance_matrix:
    input:
        pggb_input = join_path(results_dir, 'pggb', 'genomes.sample_' + str(config['sample_size']) + '_from_each_passage.fa.gz')
    output:
        split_fastas = directory(join_path(results_dir, 'split_fasta_' + str(config['sample_size']) )),
        fastani_distance_matrix = join_path(results_dir, 'fastani', 'fastani_distance_matrix.tsv'),
    params:
        list_of_files = join_path(results_dir, 'split_fasta_' + str(config['sample_size']), 'list_of_files.txt' )
    conda:
        '../envs/fastani_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        'seqkit split -O {output.split_fastas} --by-id {input.pggb_input} && '
        "find {output.split_fastas} -name '*fa.gz' -exec readlink -f {{}} \; > {params.list_of_files} && "
        'fastANI -t {threads} --fragLen 200 --ql {params.list_of_files} --rl {params.list_of_files} -o /dev/stdout  | '
        "perl -pe 's|/.*?id_||g;s|.fa.gz||g' | awk -v OFS='\\t' '{{print $1,$2,$3}}' >{output.fastani_distance_matrix}"
# Fastani:1 ends here

# [[file:../../main.org::*PGGB][PGGB:1]]
rule pggb_pangenome:
    input:
        pggb_input = join_path(results_dir, 'pggb', 'genomes.sample_' + str(config['sample_size']) + '_from_each_passage.fa.gz')
    output:
        pggb_out = directory(join_path(results_dir, 'pggb', 'sample_' + str(config['sample_size']) ))
    params:
        **config['params']['pggb']
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "n_mappings=$( zgrep -c '>' {input.pggb_input} ) && "
        "pggb -m -p {params.map_pct_id} -n $n_mappings -s {params.segment_length} -l {params.block_length} -k {params.min_match_len} -B {params.transclose_batch} -t {threads} -o {output.pggb_out} -i {input.pggb_input}"
# PGGB:1 ends here

# [[file:../../main.org::*odgi distance matrix][odgi distance matrix:1]]
rule get_distance_metrics:
    input:
        pggb_out = join_path(results_dir, 'pggb', 'sample_' + str(config['sample_size']) )
    output:
        distance_tsv = join_path(results_dir, 'pggb', 'distance_matrix.sample.' + str(config['sample_size']) + '.tsv' )
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "odgi paths -t {threads} -d -i {input.pggb_out}/*.smooth.final.og > {output.distance_tsv}"
# odgi distance matrix:1 ends here

# [[file:../../main.org::*Plot FASTANI][Plot FASTANI:1]]
rule plot_fast_ani:
    input:
        fastani_distance_matrix = join_path(results_dir, 'fastani', 'fastani_distance_matrix.tsv'),
        codes = join_path('data', 'tables', 'codes.txt'),
        script_fix_id = join_path(snakefile_path, 'scripts', 'fix_ids.py'),
        script_phylogeny_fastani = join_path(snakefile_path, 'scripts', 'phylogeny_fastani.R'),
    output:
        fastani_distance_matrixi_id_fixed = join_path(results_dir, 'fastani', 'fastani_distance_matrix.tsv'.replace('.tsv', 'ids_fixed.tsv')),
        rectangular = join_path(results_dir, 'fastani', 'ggtree.ecoli.phages.passages.rectangular.pdf'),
        daylight = join_path(results_dir, 'fastani', 'ggtree.ecoli.phages.passages.daylight.pdf'),
    conda:
        '../envs/R_env.yaml'
    shell:
        'python3 {input.script_fix_id} {input.fastani_distance_matrix} {input.codes} > {output.fastani_distance_matrixi_id_fixed} && '
        'Rscript {input.script_phylogeny_fastani} {output.fastani_distance_matrixi_id_fixed} {input.codes} {output.rectangular}'
# Plot FASTANI:1 ends here
