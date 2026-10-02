#!/usr/bin/env snakemake -s
# Snakefile for the FK2-AdFe-j13-LSK louse genome assembly pipeline


configfile: "config/config.yaml"


from pathlib import Path

if not config.get("dorado_model"):
    raise ValueError(
        "Set dorado_model in config/config.yaml or with --config dorado_model=..."
    )

POD5_DIRS = {
    f"{pod5_dir.parent.name}/{pod5_dir.parent.name}/{pod5_dir.name}": pod5_dir
    for sample_root in sorted(
        Path(config["raw_data_dir"]).glob(config["run_dir_prefix"] + "-LSK*/")
    )
    for pod5_dir in (
        sorted(sample_root.rglob("pod5"))
        + sorted(sample_root.rglob("pod5_skip"))
    )
}

if not POD5_DIRS:
    raise ValueError("No pod5 or pod5_skip directories found for FK2-AdFe-j13-LSK")

POD5_FILE_COUNT = sum(
    len(list(pod5_dir.glob("*.pod5"))) for pod5_dir in POD5_DIRS.values()
)

LOUSE = config["run_dir_prefix"]
DEVICE = config["dorado_device"]
DORADO_BIN = "resources/dorado/bin/dorado"
BUSCO_LINEAGE = config["busco_lineage"]
KMER_K = config["kmer_k"]
REF_ACCESSION = config.get("quast_reference_accession")
QUAST_REF = (
    f"results/reference/{REF_ACCESSION}/reference.fasta" if REF_ACCESSION else None
)
KRAKEN_DB = config.get("kraken2_db")


def _print_status(*lines):
    width = max(map(len, lines))
    border = f"+{'-' * (width + 2)}+"
    print("\n" + border)
    for line in lines:
        print(f"| {line:<{width}} |")
    print(border)


onstart:
    _print_status(
        "L. salmonis genome assembly",
        "WORKFLOW STARTING",
        f"Sample: {LOUSE}",
        f"POD5 signal sources: {len(POD5_DIRS)}",
        f"POD5 files discovered: {POD5_FILE_COUNT}",
    )


onerror:
    _print_status(
        "WORKFLOW FAILED",
        "Check the failed rule above and its log under logs/.",
    )


onsuccess:
    _print_status(
        "WORKFLOW FINISHED SUCCESSFULLY",
        f"Sample: {LOUSE}",
        "Results: results/",
    )


def pod5_files(wildcards):
    pod5_dir = POD5_DIRS[f"{wildcards.sample}/{wildcards.run}/{wildcards.signal}"]
    files = sorted(pod5_dir.glob("*.pod5"))
    if not files:
        raise ValueError(f"No POD5 files found in {pod5_dir}")
    return [str(file) for file in files]


rule all:
    input:
        f"results/polished/{LOUSE}/polished_assembly.fasta",
        f"results/qc/quast/{LOUSE}/report.txt",
        f"results/qc/busco/{LOUSE}/busco/short_summary.specific.{BUSCO_LINEAGE}.busco.txt",
        f"results/qc/genomescope/{LOUSE}/summary.txt",
        f"results/qc/smudgeplot/{LOUSE}/smudgeplot_smudgeplot.png",
        [f"results/qc/nanoplot/{key}/NanoPlot-report.html" for key in POD5_DIRS],
        f"results/qc/nanoplot/{LOUSE}/merged/NanoPlot-report.html",
        f"results/qc/coverage/{LOUSE}/coverage.mosdepth.summary.txt",
        *([f"results/qc/kraken2/{LOUSE}/kraken2_report.txt"] if KRAKEN_DB else []),
        *([f"results/qc/dotplot/{LOUSE}/dotplot.png"] if REF_ACCESSION else []),


rule basecall_all:
    input:
        [f"results/basecalling/{key}/basecalled.bam" for key in POD5_DIRS],


if config["dorado_install"] == "conda":

    # Wrapper pins the conda env's dorado so piped rules in other envs can call it.
    rule setup_dorado:
        output:
            DORADO_BIN,
        log:
            "logs/setup/dorado.log",
        conda:
            "workflow/envs/dorado.yaml"
        message:
            "Setting up dorado using conda environment"
        shell:
            'mkdir -p "$(dirname {output:q})" && '
            'printf \'#!/bin/sh\\nexec %s "$@"\\n\' "$(command -v dorado)" > {output:q} '
            "2> {log:q} && chmod +x {output:q}"

else:

    rule setup_dorado:
        output:
            DORADO_BIN,
        log:
            "logs/setup/dorado.log",
        params:
            version=config["dorado_version"],
        message:
            f"Setting up dorado version {config['dorado_version']}"
        shell:
            'case "$(uname -s)-$(uname -m)" in '
            "Linux-x86_64) platform=linux-x64 ;; "
            "Linux-aarch64) platform=linux-arm64 ;; "
            "Darwin-arm64) platform=osx-arm64 ;; "
            '*) echo "unsupported platform: $(uname -s)-$(uname -m)" >&2; exit 1 ;; '
            "esac && mkdir -p resources/dorado && "
            'curl -fsSL "https://cdn.oxfordnanoportal.com/software/analysis/'
            'dorado-{params.version}-$platform.tar.gz" 2> {log:q} '
            "| tar -xz -C resources/dorado --strip-components=1 2>> {log:q}"


rule basecall:
    input:
        pod5=pod5_files,
        dorado=DORADO_BIN,
    output:
        "results/basecalling/{sample}/{run}/{signal}/basecalled.bam",
    log:
        "logs/basecalling/{sample}/{run}/{signal}.log",
    params:
        model=config["dorado_model"],
        # Directing dorado to store/load models locally
        models_dir="resources/dorado_models", 
        pod5_dir=lambda wildcards: str(
            POD5_DIRS[f"{wildcards.sample}/{wildcards.run}/{wildcards.signal}"]
        ),

    shell:
        "mkdir -p {params.models_dir} && "
        "{input.dorado:q} basecaller -x 'cuda:all' --emit-moves {params.model:q} "
        "--models-directory {params.models_dir:q} {params.pod5_dir:q} > {output:q} 2> {log:q}"


# All runs of this louse are pooled into a single call set for assembly.
rule merge_calls:
    input:
        [f"results/basecalling/{key}/basecalled.bam" for key in POD5_DIRS],
    output:
        f"results/calls/{LOUSE}/calls.bam",
    log:
        f"logs/calls/{LOUSE}/merge_calls.log",
    conda:
        "workflow/envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools merge -@ {threads} -o {output:q} {input:q} 2> {log:q}"


rule calls_fastq:
    input:
        f"results/calls/{LOUSE}/calls.bam",
    output:
        f"results/calls/{LOUSE}/calls.fastq.gz",
    log:
        f"logs/calls/{LOUSE}/calls_fastq.log",
    conda:
        "workflow/envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools fastq -@ {threads} {input:q} 2> {log:q} | gzip > {output:q}"


rule correct_reads:
    input:
        reads=f"results/calls/{LOUSE}/calls.fastq.gz",
        dorado=DORADO_BIN,
    output:
        f"results/corrected/{LOUSE}/corrected.fasta",
    log:
        f"logs/corrected/{LOUSE}/correct.log",
    threads: config["threads"]
    params:
        device=DEVICE,
    shell:
        "{input.dorado:q} correct -x {params.device:q} -t {threads} {input.reads:q} "
        "> {output:q} 2> {log:q}"


rule hifiasm_assemble:
    input:
        f"results/corrected/{LOUSE}/corrected.fasta",
    output:
        f"results/assembly/{LOUSE}/hifiasm.bp.p_ctg.gfa",
    log:
        f"logs/assembly/{LOUSE}/hifiasm.log",
    conda:
        "workflow/envs/hifiasm.yaml"
    threads: config["threads"]
    params:
        prefix=f"results/assembly/{LOUSE}/hifiasm",
        extra=config["hifiasm_extra"],
    shell:
        "hifiasm -o {params.prefix:q} -t {threads} {params.extra} {input:q} "
        "> {log:q} 2>&1"


rule assembly_fasta:
    input:
        f"results/assembly/{LOUSE}/hifiasm.bp.p_ctg.gfa",
    output:
        f"results/assembly/{LOUSE}/draft_assembly.fasta",
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input:q} >{output:q}
        """


rule align_calls:
    input:
        draft=f"results/assembly/{LOUSE}/draft_assembly.fasta",
        calls=f"results/calls/{LOUSE}/calls.bam",
        dorado=DORADO_BIN,
    output:
        f"results/polished/{LOUSE}/aligned_calls.bam",
    log:
        f"logs/polished/{LOUSE}/aligner.log",
    threads: config["threads"]
    shell:
        "{input.dorado:q} aligner -t {threads} {input.draft:q} {input.calls:q} "
        "> {output:q} 2> {log:q}"


# dorado polish needs the alignments sorted and indexed by reference position.
rule sort_aligned_calls:
    input:
        f"results/polished/{LOUSE}/aligned_calls.bam",
    output:
        bam=f"results/polished/{LOUSE}/aligned_calls.sorted.bam",
        bai=f"results/polished/{LOUSE}/aligned_calls.sorted.bam.bai",
    log:
        f"logs/polished/{LOUSE}/sort_aligned_calls.log",
    conda:
        "workflow/envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} -o {output.bam:q} {input:q} 2> {log:q} && "
        "samtools index -@ {threads} {output.bam:q} 2>> {log:q}"


rule polish_assembly:
    input:
        bam=f"results/polished/{LOUSE}/aligned_calls.sorted.bam",
        bai=f"results/polished/{LOUSE}/aligned_calls.sorted.bam.bai",
        draft=f"results/assembly/{LOUSE}/draft_assembly.fasta",
        dorado=DORADO_BIN,
    output:
        f"results/polished/{LOUSE}/polished_assembly.fasta",
    log:
        f"logs/polished/{LOUSE}/polish.log",
    threads: config["threads"]
    params:
        device=DEVICE,
    shell:
        "{input.dorado:q} polish -x {params.device:q} -t {threads} {input.bam:q} "
        "{input.draft:q} > {output:q} 2> {log:q}"


rule download_reference:
    output:
        "results/reference/{accession}/reference.fasta",
    log:
        "logs/reference/{accession}/download.log",
    conda:
        "workflow/envs/ncbi-datasets.yaml"
    params:
        outdir="results/reference/{accession}",
    shell:
        "datasets download genome accession {wildcards.accession:q} "
        "--include genome --filename {params.outdir:q}/ncbi.zip > {log:q} 2>&1 && "
        "unzip -o {params.outdir:q}/ncbi.zip -d {params.outdir:q} >> {log:q} 2>&1 && "
        "cat {params.outdir:q}/ncbi_dataset/data/*/*.fna > {output:q} && "
        "rm -r {params.outdir:q}/ncbi.zip {params.outdir:q}/ncbi_dataset"


rule quast:
    input:
        assembly=f"results/polished/{LOUSE}/polished_assembly.fasta",
        reference=[QUAST_REF] if QUAST_REF else [],
    output:
        f"results/qc/quast/{LOUSE}/report.txt",
    log:
        f"logs/qc/{LOUSE}/quast.log",
    conda:
        "workflow/envs/quast.yaml"
    threads: config["threads"]
    params:
        outdir=f"results/qc/quast/{LOUSE}",
        ref_arg=f"-r {QUAST_REF}" if QUAST_REF else "",
    shell:
        "quast -o {params.outdir:q} -t {threads} {params.ref_arg} "
        "{input.assembly:q} > {log:q} 2>&1"


rule busco:
    input:
        f"results/polished/{LOUSE}/polished_assembly.fasta",
    output:
        f"results/qc/busco/{LOUSE}/busco/short_summary.specific.{BUSCO_LINEAGE}.busco.txt",
    log:
        f"logs/qc/{LOUSE}/busco.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: config["threads"]
    params:
        out_path=f"results/qc/busco/{LOUSE}",
        lineage=BUSCO_LINEAGE,
    shell:
        "busco -f -i {input:q} -m genome -l {params.lineage:q} "
        "-o busco --out_path {params.out_path:q} -c {threads} > {log:q} 2>&1"


rule genomescope_kmers:
    input:
        f"results/calls/{LOUSE}/calls.fastq.gz",
    output:
        histo=f"results/kmers/{LOUSE}/genomescope.histo",
        jf=temp(f"results/kmers/{LOUSE}/genomescope.jf"),
    log:
        f"logs/kmers/{LOUSE}/genomescope_kmers.log",
    conda:
        "workflow/envs/genomescope.yaml"
    threads: config["threads"]
    params:
        k=KMER_K,
    shell:
        "zcat {input:q} | jellyfish count -C -m {params.k} -s 1G -t {threads} "
        "-o {output.jf:q} /dev/stdin 2> {log:q} && "
        "jellyfish histo -t {threads} {output.jf:q} > {output.histo:q} 2>> {log:q}"


rule genomescope:
    input:
        f"results/kmers/{LOUSE}/genomescope.histo",
    output:
        f"results/qc/genomescope/{LOUSE}/summary.txt",
    log:
        f"logs/qc/{LOUSE}/genomescope.log",
    conda:
        "workflow/envs/genomescope.yaml"
    params:
        outdir=f"results/qc/genomescope/{LOUSE}",
        k=KMER_K,
        ploidy=config["genomescope_ploidy"],
    shell:
        "genomescope2 -i {input:q} -o {params.outdir:q} -k {params.k} "
        "-p {params.ploidy} > {log:q} 2>&1"


rule smudgeplot_kmers:
    input:
        f"results/calls/{LOUSE}/calls.fastq.gz",
    output:
        pre=f"results/kmers/{LOUSE}/kmcdb.kmc_pre",
        suf=f"results/kmers/{LOUSE}/kmcdb.kmc_suf",
        hist=f"results/kmers/{LOUSE}/kmcdb.hist",
    log:
        f"logs/kmers/{LOUSE}/smudgeplot_kmers.log",
    conda:
        "workflow/envs/smudgeplot.yaml"
    threads: config["threads"]
    params:
        k=KMER_K,
        mem=config["kmc_memory_gb"],
        prefix=f"results/kmers/{LOUSE}/kmcdb",
        tmpdir=f"results/kmers/{LOUSE}/kmc_tmp",
    shell:
        "mkdir -p {params.tmpdir:q} && "
        "kmc -k{params.k} -t{threads} -m{params.mem} -ci1 -cs10000 "
        "{input:q} {params.prefix:q} {params.tmpdir:q} > {log:q} 2>&1 && "
        "kmc_tools transform {params.prefix:q} histogram {output.hist:q} "
        "-cx10000 >> {log:q} 2>&1 && "
        "rm -r {params.tmpdir:q}"


rule smudgeplot:
    input:
        pre=f"results/kmers/{LOUSE}/kmcdb.kmc_pre",
        suf=f"results/kmers/{LOUSE}/kmcdb.kmc_suf",
        hist=f"results/kmers/{LOUSE}/kmcdb.hist",
    output:
        f"results/qc/smudgeplot/{LOUSE}/smudgeplot_smudgeplot.png",
    log:
        f"logs/qc/{LOUSE}/smudgeplot.log",
    conda:
        "workflow/envs/smudgeplot.yaml"
    params:
        db=f"results/kmers/{LOUSE}/kmcdb",
        dump=f"results/qc/smudgeplot/{LOUSE}/kmcdb.dump",
        pair_prefix=f"results/qc/smudgeplot/{LOUSE}/kmer_pairs",
        out_prefix=f"results/qc/smudgeplot/{LOUSE}/smudgeplot",
    shell:
        "L=$(smudgeplot.py cutoff {input.hist:q} L) && "
        "U=$(smudgeplot.py cutoff {input.hist:q} U) && "
        'kmc_tools transform {params.db:q} -ci"$L" -cx"$U" dump -s {params.dump:q} '
        "> {log:q} 2>&1 && "
        "smudgeplot.py hetkmers -o {params.pair_prefix:q} < {params.dump:q} "
        ">> {log:q} 2>&1 && "
        "smudgeplot.py plot {params.pair_prefix:q}_coverages.tsv -o {params.out_prefix:q} "
        ">> {log:q} 2>&1 && "
        "rm {params.dump:q}"


# Dotplot of the polished assembly against the configured NCBI reference.
rule dotplot:
    input:
        assembly=f"results/polished/{LOUSE}/polished_assembly.fasta",
        reference=[QUAST_REF] if QUAST_REF else [],
    output:
        f"results/qc/dotplot/{LOUSE}/dotplot.png",
    log:
        f"logs/qc/{LOUSE}/dotplot.log",
    conda:
        "workflow/envs/mummer.yaml"
    threads: config["threads"]
    params:
        prefix=f"results/qc/dotplot/{LOUSE}/dotplot",
    shell:
        "nucmer -t {threads} -p {params.prefix:q} {input.reference:q} "
        "{input.assembly:q} > {log:q} 2>&1 && "
        "delta-filter -1 {params.prefix:q}.delta > {params.prefix:q}.1delta "
        "2>> {log:q} && "
        "mummerplot --fat --png -p {params.prefix:q} {params.prefix:q}.1delta "
        ">> {log:q} 2>&1"


# Per-run read QC to spot bad flowcells or libraries before pooling.
rule nanoplot_run:
    input:
        "results/basecalling/{sample}/{run}/{signal}/basecalled.bam",
    output:
        "results/qc/nanoplot/{sample}/{run}/{signal}/NanoPlot-report.html",
    log:
        "logs/qc/nanoplot/{sample}/{run}/{signal}.log",
    conda:
        "workflow/envs/nanoplot.yaml"
    threads: config["threads"]
    params:
        outdir="results/qc/nanoplot/{sample}/{run}/{signal}",
    shell:
        "NanoPlot -t {threads} --ubam {input:q} -o {params.outdir:q} " "> {log:q} 2>&1"


rule nanoplot_merged:
    input:
        f"results/calls/{LOUSE}/calls.fastq.gz",
    output:
        f"results/qc/nanoplot/{LOUSE}/merged/NanoPlot-report.html",
    log:
        f"logs/qc/nanoplot/{LOUSE}/merged.log",
    conda:
        "workflow/envs/nanoplot.yaml"
    threads: config["threads"]
    params:
        outdir=f"results/qc/nanoplot/{LOUSE}/merged",
    shell:
        "NanoPlot -t {threads} --fastq {input:q} -o {params.outdir:q} " "> {log:q} 2>&1"


# Per-read output is discarded; the taxonomic report is what we inspect.
rule kraken2:
    input:
        f"results/calls/{LOUSE}/calls.fastq.gz",
    output:
        f"results/qc/kraken2/{LOUSE}/kraken2_report.txt",
    log:
        f"logs/qc/{LOUSE}/kraken2.log",
    conda:
        "workflow/envs/kraken2.yaml"
    threads: config["threads"]
    params:
        db=KRAKEN_DB or "",
    shell:
        "kraken2 --db {params.db:q} --threads {threads} --gzip-compressed "
        "--report {output:q} --output /dev/null {input:q} > {log:q} 2>&1"


# Coverage is assessed on the polished assembly, so calls are re-aligned to it.
rule align_calls_polished:
    input:
        assembly=f"results/polished/{LOUSE}/polished_assembly.fasta",
        calls=f"results/calls/{LOUSE}/calls.bam",
        dorado=DORADO_BIN,
    output:
        bam=f"results/qc/coverage/{LOUSE}/calls_on_polished.sorted.bam",
        bai=f"results/qc/coverage/{LOUSE}/calls_on_polished.sorted.bam.bai",
    log:
        f"logs/qc/{LOUSE}/align_calls_polished.log",
    conda:
        "workflow/envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "{input.dorado:q} aligner -t {threads} {input.assembly:q} {input.calls:q} 2> {log:q} "
        "| samtools sort -@ {threads} -o {output.bam:q} - 2>> {log:q} && "
        "samtools index -@ {threads} {output.bam:q} 2>> {log:q}"


rule coverage:
    input:
        bam=f"results/qc/coverage/{LOUSE}/calls_on_polished.sorted.bam",
        bai=f"results/qc/coverage/{LOUSE}/calls_on_polished.sorted.bam.bai",
    output:
        f"results/qc/coverage/{LOUSE}/coverage.mosdepth.summary.txt",
    log:
        f"logs/qc/{LOUSE}/mosdepth.log",
    conda:
        "workflow/envs/mosdepth.yaml"
    threads: config["threads"]
    params:
        prefix=f"results/qc/coverage/{LOUSE}/coverage",
    shell:
        "mosdepth -t {threads} -x -n --by 10000 {params.prefix:q} {input.bam:q} "
        "> {log:q} 2>&1"
