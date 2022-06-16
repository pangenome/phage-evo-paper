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

# [[file:../../main.org::*Extract phage genomes][Extract phage genomes:1]]
rule filter_by_lenght_to_get_phage_genomes:
    input:
        reads = join_path(config['data']['reads'], '{sample}' + config['data']['reads_suffix']),
    output:
        putative_phage_genomes = join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.fasta'),
    params:
        **config['params']['seqkit'],
    threads:
        get_cores_perc(0.3)
    conda:
        '../envs/graphaligner_env.yaml'
    shell:
        'seqkit seq {input.reads} -j {threads} -m {params.min} -M {params.max} | seqkit fq2fa | '
        "sed -r '/>/ s|>|>{wildcards.sample}#1#|' | bgzip > {output.putative_phage_genomes} "
# Extract phage genomes:1 ends here

# [[file:../../main.org::*Filter genomes][Filter genomes:1]]
rule filter_out_bacterial_genomes:
    input:
        target = config['data']['genomes']['ecoli_and_phages'],
        putative_phage_genomes = expand(join_path(results_dir, 'minia', '{sample}', '{sample}' + '.putative_phage_genomes' + '.fasta'), sample=SAMPLES),
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
        "cat {input.putative_phage_genomes} | bgzip -@ {threads} >{output.all_genomes_merged} && "
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
        sample_ids = join_path(results_dir, 'pggb', '{experiment}', '{experiment}.ids.sample_size_' + str(config['sample_size']) + '.txt'),
    params:
        sample_size = config['sample_size'],
        log_dir = join_path(str(Path('results').parent.absolute()), 'logs'),
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_{wildcards.experiment}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        "awk -F$'\\t' '/^{wildcards.experiment}/ {{print $3}}' {input.codes}  | "
        'while read f; do grep -P "^${{f}}#" {input.ids_to_keep} | shuf -n {params.sample_size}; done | tee {output.sample_ids} && '
        "samtools faidx {input.all_genomes_merged_filtered} -r {output.sample_ids} | "
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
        'gunzip {output.split_fastas}/*.fa.gz && '
        "find {output.split_fastas} -name '*.fa' -exec readlink -f {{}} \; > {output.list_of_files} && "
        'fastANI -t {threads} --fragLen {params.frag_lenght} --ql {output.list_of_files} --rl {output.list_of_files} -o /dev/stdout  | '
        "perl -pe 's|/.*?id_||g;s|.fa||g' | awk -v OFS='\\t' '{{print $1,$2,$3}}' >{output.fastani_distance_matrix}"
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

# [[file:../../main.org::*ORF prediction PHANOTATE][ORF prediction PHANOTATE:1]]
rule orf_prediction_phanotate:
    input:
        list_of_files = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.list_of_splited_fastas_pahts.sample_size_' + str(config['sample_size']) + '.txt'),
        phanotate_runner = join_path(scripts_dir, 'phanotate_runner.py'),
    output:
        phanotate_dir = directory(join_path(results_dir, 'annotations', 'phanotate', '{experiment}')),
        finished = join_path(results_dir, 'annotations', 'phanotate', '{experiment}', 'finished_phanotate'),
    params:
        **config['params']['phanotate'],
        log_dir = join_path(snakefile_path, '..', 'logs'),
    conda:
        '../envs/phanotate_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_{wildcards.experiment}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        'python3 {input.phanotate_runner} --input_file_list {input.list_of_files} '
        ' --threads {threads} --out_format {params.out_format} --output_dir {output.phanotate_dir} && '
        '>{output.finished} '
# ORF prediction PHANOTATE:1 ends here

# [[file:../../main.org::*PROKKA annotation][PROKKA annotation:1]]
rule annotation_prokka:
    input:
        list_of_files = join_path(results_dir, 'fastani', '{experiment}', '{experiment}.list_of_splited_fastas_pahts.sample_size_' + str(config['sample_size']) + '.txt'),
        mmseqs_phrogs_db = join_path('data', 'phrogs_db', 'phrogs_rep_seq.fasta'),
        prokka_runner = join_path(scripts_dir, 'prokka_runner.py'),
    output:
        prokka_dir = directory(join_path(results_dir, 'annotations', 'prokka', '{experiment}')),
        finished = join_path(results_dir, 'annotations', 'prokka', '{experiment}', 'finished_prokka'),
    params:
        log_dir = join_path(snakefile_path, '..', 'logs'),
    conda:
        '../envs/prokka_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_{wildcards.experiment}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        'python3 {input.prokka_runner} --input_file_list {input.list_of_files} '
        '--output_dir {output.prokka_dir} --proteins {input.mmseqs_phrogs_db} '
        '--threads {threads} --prokka_threads 4 && '
        'touch {output.finished}'
# PROKKA annotation:1 ends here

# [[file:../../main.org::*Download prokka proteins db (Phrogs)][Download prokka proteins db (Phrogs):1]]
rule download_phrogs_database:
    output:
        phrogs_tar = join_path('data', 'phrogs_db', 'FAA_phrog.tar.gz'),
    params:
        log_dir = join_path(snakefile_path, '..', 'logs'),
        mmseqs_phrogs_url = config['mmseqs_phrogs_url'],
    threads:
        1
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        'wget -O {output.phrogs_tar} {params.mmseqs_phrogs_url}'
# Download prokka proteins db (Phrogs):1 ends here

# [[file:../../main.org::*Cluster prokka proteins db (Phrogs)][Cluster prokka proteins db (Phrogs):1]]
rule reclust_phrogs_database:
    input:
        phrogs_tar = join_path('data', 'phrogs_db', 'FAA_phrog.tar.gz'),
    output:
        mmseqs_phrogs_db = join_path('data', 'phrogs_db', 'phrogs_rep_seq.fasta'),
    params:
        log_dir = join_path(snakefile_path, '..', 'logs'),
        mmseqs_phrogs_url = config['mmseqs_phrogs_url'],
        mmseqs_multifasta_dir = join_path('data', 'phrogs_db', 'FAA_phrog'),
        phrogs_db_dir = join_path('data', 'phrogs_db')
    threads:
        get_cores_perc(1)
    conda:
        '../envs/mmseqs2_env.yaml'
    shell:
        'exec &> >( tee {params.log_dir}/{rule}_$(date +%Y_%m_%d_-_%H_%M_%S).log ) && '
        'tar -xf {input.phrogs_tar} -C {params.phrogs_db_dir} && '
        'cat {params.mmseqs_multifasta_dir}/*.faa > {params.phrogs_db_dir}/multifasta.faa && '
        'mmseqs easy-cluster {params.phrogs_db_dir}/multifasta.faa {params.phrogs_db_dir}/phrogs {params.phrogs_db_dir}/tmp --threads {threads}'
# Cluster prokka proteins db (Phrogs):1 ends here
