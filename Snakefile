configfile: "config.yaml"

if config["mode"] not in ["mapper", "assembler"]:
    raise ValueError("config.yaml error: 'mode' must be either 'mapper' or 'assembler'")

rule all:
    input:
        config["output_mapping"] if config["mode"] == "mapper"
        else config["output_contigs"]

rule run_mapper:
    input:
        ref=config["input_reference"],
        reads=config["input_reads"]
    output:
        config["output_mapping"]
    conda:
        "env.yml"
    shell:
        "python3 mapper/mapper.py {input.ref} {input.reads} {output}"

rule run_assembly:
    input:
        config["input_reads"]
    output:
        config["output_contigs"]
    conda:
        "env.yml"
    shell:
        "python3 assembler/assembly {input} {output}"
