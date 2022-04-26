# [[file:../../main.org::*Merge reads][Merge reads:1]]
# rule merged_reads:
#     input:
#         expand(config['reads']+'/{sample}.fastq', sample=SAMPLES)
#     output:
#         config['data']+'/merged.before.filtering.fastq'
#     shell:
#         "cat {input} > {output}"
# Merge reads:1 ends here

# [[file:../../main.org::*Map reads to minia assembly][Map reads to minia assembly:1]]
rule prefix_fastq:
    input:
        samples=expand(join_path(config['data']['reads'], '{sample}.merged.fastq'), sample=SAMPLES),
    params:
        samples_prefixed=join_path(config['data']['reads'], 'P1-10.merged.prefixed.fastq')
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
        cat {params.samples_prefixed} | pigz -p {threads} > {output.samples_prefixed_gzipped}
        """
# Map reads to minia assembly:1 ends here

# [[file:../../main.org::*NANOPLOT][NANOPLOT:1]]
rule nanoplot:
    input:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.{state}_qc.fastq.gz'),
    output:
        directory("outputs/nanoplot/{state}_filter")
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
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.before_qc.fastq.gz'),
    output:
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.after_qc.fastq.gz'),
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
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.after_qc.fastq.gz'),
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
        samples_prefixed_gzipped=join_path(config['data']['reads'], 'P1-10.merged.prefixed.before_qc.fastq.gz'),
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
        minia_assembly_polished_filtered = filter_contigs_prefix + '.contigs.polished.fa'
    params:
        **config['params']['minia']
    conda:
        'envs/bio_env.yaml'
    shell:
        "python3 {input.script} {params.min_contig_lenght}  {params.max_contig_lenght} > {output.minia_assembly_polished_filtered}"
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
        'envs/pggb_env.yaml'
    shell:
        "cat {input.minia_assembly_polished_filtered} | bzip -@ {threads} > {output.minia_assembly_polished_filtered_crompressed}"
        "samtools faidx {output.minia_assembly_polished_filtered_crompressed}"
# Create index:1 ends here

# [[file:../../main.org::*Get sample and add parental phages genomes][Get sample and add parental phages genomes:1]]
# rule add_parental_genomes_and_get_sample:
#     input:
# Get sample and add parental phages genomes:1 ends here

# [[file:../../main.org::*PGGB minia_polished][PGGB minia_polished:1]]
rule pggb_minia:
    input:
        minia_assembly_polished_filtered_crompressed = filter_contigs_prefix + '.contigs.polished.fa.gz',
        fai = filter_contigs_prefix + '.contigs.polished.fa.gz.fai',
        gzi = filter_contigs_prefix + '.contigs.polished.fa.gz.gzi',
    output:
        directory("outputs/pggb/minia"+pggb_prefix),
    params:
        **config['params']['pggb']
    conda:
        'envs/pggb_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        "n_mappings=$( grep -c '>' {input.minia_assembly_polished_filtered_crompressed} )" # get number of genomes
        "pggb -m -p {params.map_pct_id} -n $n_mappings -s {params.segment_length} -l {params.block_lengt} -t {threads} -o {output} -i {input.minia_assembly_polished_filtered}"
# PGGB minia_polished:1 ends here