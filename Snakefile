rule check_setup:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/check_setup.R"
rule filter_reads:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/filtering_reads.R"

rule estimate_errors:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/estimate_errors.R"

rule infer_sequences:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/infer_sequences.R"

rule merge_sequences:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/merge_sequences.R"

rule make_sequence_table:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/make_sequence_table.R"

rule remove_chimeras:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/remove_chimeras.R"

rule assign_taxonomy:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/assign_taxonomy.R"

rule align_sequences:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/align_sequences.R"

rule make_tree:
    input:
        "data/output/run_1/temp/alignment.fasta"
    output:
        "data/output/run_1/temp/fasttree.tree"
    shell:
        "FastTree -nt {input} > {output}"

rule make_phyloseq:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/make_phyloseq.R"

rule track_reads:
    input:
        config_file="config.csv",
        extra_functions="scripts/functions.R"
    script:
        "scripts/track_reads.R"