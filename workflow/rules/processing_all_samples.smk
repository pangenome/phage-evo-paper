# [[file:../../main.org::*Plot lengths][Plot lengths:1]]
rule quality_check_plot:
    input:
        reads = join_path(config['data']['reads'], '{sample}' + config['data']['reads_suffix'])
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
        reads = join_path(config['data']['reads'], '{sample}' + config['data']['reads_suffix']),
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
        reads = join_path(config['data']['reads'], '{sample}' + config['data']['reads_suffix']),
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

# [[file:../../main.org::*Filter genomes][Filter genomes:1]]
rule filter_out_bacterial_genomes:
    input:
        target = config['data']['genomes']['ecoli_and_phages'],
        putative_phage_genomes_polished = expand(join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.polished' + '.prefixed' + '.fa.gz'), sample=SAMPLES),
    output:
        all_genomes_merged = join_path(results_dir, 'pggb', 'all_genomes_merged.fa.gz'),
        all_genomes_merged_filtered = join_path(results_dir, 'pggb', 'all_genomes_merged.filter_out_bacteria.fa.gz'),
        ids_to_keep = join_path(results_dir, 'pggb', 'ids_to_keep.txt'),
    params:
        **config['params']['removing_bacteria'],
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "zcat {input.putative_phage_genomes_polished} | bgzip -@ {threads} >{output.all_genomes_merged} && "
        "samtools faidx {output.all_genomes_merged} && "
        "samtools faidx {output.all_genomes_merged} "
        "-r <(wfmash {input.target} {output.all_genomes_merged} -s {params.segment_length} -l {params.block_length} -p {params.map_pct_id} -t {threads} | "
        "awk -v min_qcov={params.min_qcov} '/E_coli/ {{ qcov=$11/$2; if ( !(qcov >= min_qcov) ) print $1; }}' | sort -u | tee {output.ids_to_keep} ) > "
        "{output.all_genomes_merged_filtered}"
# Filter genomes:1 ends here

# [[file:../../main.org::*Sample genomes][Sample genomes:1]]
rule sample_genomes:
    input:
        all_genomes_merged_filtered = join_path(results_dir, 'pggb', 'all_genomes_merged.filter_out_bacteria.fa.gz'),
        ids_to_keep = join_path(results_dir, 'pggb', 'ids_to_keep.txt'),
        codes = join_path('data', 'tables', 'codes.tsv'),
    output:
        pggb_input = join_path(results_dir, 'pggb', '{experiment}', '{experiment}.merged_genomes.sample_size_' + str(config['sample_size']) + '.fa.gz'),
    params:
        sample_size = config['sample_size'],
        log_dir = join_path(str(Path('results').parent.absolute()), 'logs'),
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_{wildcards.experiment}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        "samtools faidx {input.all_genomes_merged_filtered} -r "
        "<( awk -F$'\\t' '/^{wildcards.experiment}/ {{print $3}}' {{input.codes}}  | "
        'while read f; do grep -P "^${{f}}#" {input.ids_to_keep} | shuf -n {params.sample_size}; done ) | '
        'bgzip -@ {threads} > {output.pggb_input} '
# Sample genomes:1 ends here

# [[file:../../main.org::*Fastani][Fastani:1]]
rule fastaANI_distance_matrix:
    input:
        pggb_input = join_path(results_dir, 'pggb', '{experiment}', '{experiment}.merged_genomes.sample_size_' + str(config['sample_size']) + '.fa.gz'),
    output:
        split_fastas = directory(join_path(results_dir, 'fastani',  '{experiment}', '{experiment}.split_fasta_' + str(config['sample_size']) )),
        fastani_distance_matrix = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.fastani_distance_matrix.sample_size_' + str(config['sample_size']) + '.tsv'),
        list_of_files = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.list_of_splited_fastas_pahts.sample_size_' + str(config['sample_size']) + '.txt'),
    params:
        **config['params']['fastani'],
        log_dir = join_path(str(Path('results').parent.absolute()), 'logs')
    conda:
        '../envs/fastani_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_{wildcards.experiment}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        'seqkit split -O {output.split_fastas} --by-id {input.pggb_input} && '
        "find {output.split_fastas} -name '*fa.gz' -exec readlink -f {{}} \; > {output.list_of_files} && "
        'fastANI -t {threads} --fragLen {params.frag_lenght} --ql {output.list_of_files} --rl {output.list_of_files} -o /dev/stdout  | '
        "perl -pe 's|/.*?id_||g;s|.fa.gz||g' | awk -v OFS='\\t' '{{print $1,$2,$3}}' >{output.fastani_distance_matrix}"
# Fastani:1 ends here

# [[file:../../main.org::*Plot FASTANI][Plot FASTANI:1]]
rule plot_fast_ani:
    input:
        fastani_distance_matrix = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.fastani_distance_matrix.sample_size_' + str(config['sample_size']) + '.tsv'),
        codes = join_path('data', 'tables', 'codes.tsv'),
        script_fix_id = join_path(snakefile_path, 'scripts', 'fix_ids.py'),
        script_phylogeny_fastani = join_path(snakefile_path, 'scripts', 'phylogeny_fastani.R'),
    output:
        fastani_distance_matrix_id_fixed = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.fastani_distance_matrix.sample_size_' + str(config['sample_size']) + '.ids_fixed.tsv'),
        rectangular = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.ggtree.ecoli.phages.passages.rectangular.sample_size_' + str(config['sample_size']) +  '.pdf'),
    params:
        title = "{}.sample_size_{}.K{}.A{}_bp_relative.min40K.max50.GA_polished".format('{experiment}', config['sample_size'], config['params']['minia']['kmer'], config['params']['minia']['P1_abundance'])
    conda:
        '../envs/R_env.yaml'
    shell:
        'python3 {input.script_fix_id} {input.fastani_distance_matrix} {input.codes} > {output.fastani_distance_matrix_id_fixed} && '
        'Rscript {input.script_phylogeny_fastani} {output.fastani_distance_matrix_id_fixed} {input.codes} {output.rectangular} {params.title}'
# Plot FASTANI:1 ends here
