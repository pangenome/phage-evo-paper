# [[file:../../main.org::*Map reads to minia assembly][Map reads to minia assembly:1]]
rule prefix_fastq:
    input:
        sample=join_path(config['data']['reads'], '{sample}.merged.fastq'),
    output:
        sample_prefixed=join_path(config['data']['reads'], '{sample}.merged.prefixed.fastq')
    threads:
        1
    shell:
        "sed -r '/^@/ s/^@/@{wildcards.sample}#1#/' {input.sample} > {output.sample_prefixed}"
# Map reads to minia assembly:1 ends here
