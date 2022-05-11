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
