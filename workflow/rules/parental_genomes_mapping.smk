# [[file:../../main.org::*Map reads to minia assembly][Map reads to minia assembly:1]]
rule prefix_fastq:
    input:
        samples=expand(join_path(config['data']['reads'], '{sample}.merged.fastq'), sample=SAMPLES),
    output:
        samples_prefixed=join_path(config['data']['reads'], 'P1-10.merged.prefixed.fastq')
    threads:
        1
    shell:
        """
        echo {input.samples} \
            | tr ' ' '\\n' \
            | while read sample; do
                prefix=$( basename $sample | cut -d'.' -f1)
                sed -r '/^@.+runid/ s/^@/@'$prefix'#1#/' $sample >> {output.samples_prefixed}
            done
        """
# Map reads to minia assembly:1 ends here
