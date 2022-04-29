# [[file:../../main.org::*Prefix reads][Prefix reads:1]]
rule prefix_fastq:
    input:
        sample = join_path(config['data']['reads'], '{sample}.merged.fastq'),
    output:
        sample_prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz')
    threads:
        get_cores_perc(1)
    conda:
        '../envs/pggb_env.yaml'
    shell:
        "prefix=$( basename {input.sample} | cut -d'.' -f1) && "
        "sed -r '/^@.+runid/ s/^@/@'$prefix'#1#/' {input.sample} | bgzip > {output.sample_prefixed}"
# Prefix reads:1 ends here

# [[file:../../main.org::*nanoplot][nanoplot:1]]
rule quality_check_plot_before_filtering:
    input:
        sample = join_path(config['data']['reads'], '{state}', '{sample}.{state}.fastq.gz')
    output:
        plot_dir = directory("results/single/nanoplot/{state}/{sample}")
    threads:
        get_cores_perc(1)
    conda:
        "../envs/nanoplot_env.yaml"
    shell:
        "NanoPlot -t 2 --plots dot -o {output.plot_dir} --fastq {input.sample}"
# nanoplot:1 ends here

# [[file:../../main.org::*FILTER READS][FILTER READS:1]]
rule filter_reads:
    input:
        prefixed = join_path(config['data']['reads'], 'prefixed', '{sample}.prefixed.fastq.gz')
    output:
        filtered = join_path(config['data']['reads'], 'filtered', '{sample}.filtered.fastq.gz')
    params:
        **config['params']['filtlong']
    conda:
        "../envs/filtlong_env.yaml"
    threads:
        get_cores_perc(1)
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input.prefixed} | pigz > {output.filtered}"
# FILTER READS:1 ends here
