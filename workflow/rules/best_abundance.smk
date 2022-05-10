# [[file:../../main.org::*Minia assembly][Minia assembly:1]]
rule minia:
    input:
        filtered = join_path(config['data']['reads'], 'prefixed', 'P1.prefixed.fastq.gz')
    output:
        minia_assembly =  join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly.contigs.fa')
    threads:
        6
    params:
        kmer = '{kmer}',
        abundance = '{abundance}',
    conda:
        '../envs/minia_env.yaml'
    shell:
        "minia -nb-cores 5 -kmer-size {params.kmer} -abundance-min {params.abundance} -out $( echo {output.minia_assembly} | sed 's/.contigs.fa//' ) -in {input.filtered} && "
        " find $( dirname {output.minia_assembly} ) -type f ! -name '*'$(basename {output.minia_assembly}) -exec rm {{}} \;"
# Minia assembly:1 ends here

# [[file:../../main.org::*fasta to gfa][fasta to gfa:1]]
rule minia_fasta_to_gfa:
    input:
        minia_assembly =  join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly.contigs.fa'),
        script = join_path(snakefile_path, 'scripts', 'convertToGFA.py'),
    output:
        minia_assembly_gfa =  join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly.contigs.gfa'),
    params:
        kmer = '{kmer}',
    conda:
        '../envs/minia_env.yaml'
    threads:
        10
    shell:
        "python {input.script} {input.minia_assembly} {output.minia_assembly_gfa} {params.kmer}"
# fasta to gfa:1 ends here

# [[file:../../main.org::*Graphaligner MINIA][Graphaligner MINIA:1]]
# rule polishing_graphaligner_minia:
#     input:
#         samples_prefixed_gzipped = join_path(config['data']['reads'], 'prefixed', 'P1.prefixed.fastq.gz'),
#         minia_assembly_gfa =  join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly', '.contigs.gfa'),
#     output:
#         minia_gaf = join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly', '.contigs.gaf'),
#         minia_assembly_gfa_polished = join_path('results', 'test_abundance', 'minia', 'A{abundance}.K{kmer}', 'minia.assembly', '.contigs.polished.fa'),
#     threads:
#         4
#     params:
#         dbtype = "vg",
#         seed_minimizer = 15
#     conda:
#         '../envs/graphaligner_env.yaml'
#     shell:
#         "GraphAligner -g {input.minia_assembly_gfa} -f {input.samples_prefixed_gzipped} -x {params.dbtype} --threads 10 --seeds-minimizer-length {params.seed_minimizer} --seeds-minimizer-windowsize {params.seed_minimizer} -a {output.minia_gaf} --corrected-out {output.minia_assembly_gfa_polished}"
# Graphaligner MINIA:1 ends here
