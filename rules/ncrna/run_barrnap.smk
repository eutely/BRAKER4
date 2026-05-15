"""
rRNA gene prediction (via pybarrnap) and merging with BRAKER annotations.

Runs pybarrnap on the genome to predict ribosomal RNA genes (18S, 28S,
5.8S, 5S), then merges the rRNA features (along with tRNA, Infernal, and
lncRNA results when also enabled) into the final GFF3 annotation.

pybarrnap (https://github.com/moshi4/pybarrnap) is a Python re-implementation
of Seemann's barrnap. It uses Rfam 14.10 HMM profiles via pyhmmer (no
external nhmmer dependency) and resolves the eukaryotic 5.8S/28S overlap
that barrnap v0.9 produces (tseemann/barrnap#47).

Only included when run_ncrna = 1 in config.ini.

Container: quay.io/biocontainers/pybarrnap:0.5.1--pyhdfd78af_0
"""


rule run_barrnap:
    """Run pybarrnap to predict rRNA genes on the genome."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        gff="output/{sample}/barrnap/rRNA.gff3"
    log:
        "logs/{sample}/barrnap/barrnap.log"
    benchmark:
        "benchmarks/{sample}/barrnap/barrnap.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BARRNAP_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff})

        pybarrnap --kingdom euk --threads {threads} \
            {input.genome} > {output.gff}.tmp 2> {log} || true

        # Prefix rRNA gene IDs with sample name and add unique counter
        if [ -s {output.gff}.tmp ] && grep -qv '^#' {output.gff}.tmp; then
            # POSIX-portable Name extraction: the biocontainer ships busybox
            # awk, which does NOT support gawk's 3-arg match($str, regex, array).
            # Use 2-arg match() with RSTART/RLENGTH + substr() to capture the
            # old Name value instead.
            awk -F'\t' -v OFS='\t' -v p="{wildcards.sample}" '
                BEGIN {{n=1}}
                /^#/ {{print; next}}
                {{
                    oldname = ""
                    if (match($9, /Name=[^;]+/)) {{
                        # Skip the literal "Name=" prefix (5 chars).
                        oldname = substr($9, RSTART + 5, RLENGTH - 5)
                    }}
                    newname = p "-rRNA_" n "_" oldname
                    gsub(/Name=[^;]+/, "Name=" newname, $9)
                    if ($9 ~ /ID=/) {{
                        gsub(/ID=[^;]+/, "ID=" newname, $9)
                    }} else {{
                        $9 = "ID=" newname ";" $9
                    }}
                    n++
                    print
                }}
            ' {output.gff}.tmp > {output.gff}
        else
            echo "##gff-version 3" > {output.gff}
            echo "No rRNA genes found" >> {log}
        fi
        rm -f {output.gff}.tmp

        n_rrna=$(grep -cv '^#' {output.gff} || echo 0)
        echo "pybarrnap predicted $n_rrna rRNA features" >> {log}

        # Record software version (LC_ALL=C avoids locale warnings in biocontainer)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        PYBARRNAP_VER=$(LC_ALL=C pybarrnap --version 2>&1 | sed -n 's/.*\([0-9][0-9]*\.[0-9][0-9]*[^ ]*\).*/\1/p' | head -1 || true)
        ( flock 9; printf "pybarrnap\t%s\n" "$PYBARRNAP_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite pybarrnap "$REPORT_DIR" || true
        """


def _get_merge_ncrna_inputs(wildcards):
    """Get inputs for merging rRNA, tRNA, Infernal, and lncRNA into GFF3."""
    sample = wildcards.sample
    inputs = {
        "gff3": f"output/{sample}/braker.gff3",
        "rrna": f"output/{sample}/barrnap/rRNA.gff3",
        "trna": f"output/{sample}/ncrna/tRNAs.gff3",
        "infernal": f"output/{sample}/ncrna/ncRNAs_infernal.gff3",
    }
    # lncRNA only when transcript evidence is available for this sample
    mode = get_braker_mode(sample)
    if mode in ('et', 'etp', 'isoseq', 'dual'):
        inputs["lncrna"] = f"output/{sample}/ncrna/lncRNAs.gff3"
    return inputs


rule merge_rrna_into_gff3:
    """Merge barrnap rRNA, tRNA, Infernal, and lncRNA into the final GFF3."""
    input:
        unpack(_get_merge_ncrna_inputs)
    output:
        merged="output/{sample}/braker_with_ncRNA.gff3"
    log:
        "logs/{sample}/barrnap/merge_rrna.log"
    benchmark:
        "benchmarks/{sample}/merge_rrna_into_gff3/merge_rrna_into_gff3.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        """
        set -euo pipefail
        # Copy the BRAKER GFF3 as base
        cp {input.gff3} {output.merged}

        # Append rRNA features (skip header lines from barrnap output)
        grep -v '^#' {input.rrna} >> {output.merged} 2>/dev/null || true
        n_rrna=$(grep -cv '^#' {input.rrna} || echo 0)
        echo "Appended $n_rrna rRNA features" > {log}

        # Append ncRNA features if run_ncrna is enabled
        for ncrna_file in output/{wildcards.sample}/ncrna/tRNAs.gff3 \
                          output/{wildcards.sample}/ncrna/ncRNAs_infernal.gff3 \
                          output/{wildcards.sample}/ncrna/lncRNAs.gff3; do
            if [ -f "$ncrna_file" ]; then
                n_before=$(grep -cv '^#' {output.merged} || echo 0)
                grep -v '^#' "$ncrna_file" >> {output.merged} 2>/dev/null || true
                n_after=$(grep -cv '^#' {output.merged} || echo 0)
                n_added=$((n_after - n_before))
                echo "Appended $n_added features from $(basename $ncrna_file)" >> {log}
            fi
        done

        n_braker=$(grep -cv '^#' {input.gff3} || echo 0)
        n_total=$(grep -cv '^#' {output.merged} || echo 0)
        echo "Total: $n_braker BRAKER + $((n_total - n_braker)) ncRNA = $n_total features" >> {log}
        """
