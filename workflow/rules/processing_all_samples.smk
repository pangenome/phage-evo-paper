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
