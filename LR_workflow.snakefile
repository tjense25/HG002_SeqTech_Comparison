from os.path import join
import pandas as pd
import glob
tmpdir= config["scratch_dir"]
outdir=config["output_dir"]

sample_table= pd.read_table(config["sample_table"])

fastq_dict={row.SAMPLE : row.FASTQ_PATH for i,row in sample_table.iterrows()}
model_map={row.SAMPLE : row.MODEL for i,row in sample_table.iterrows()}

samples = list(fastq_dict.keys())


localrules: all
               
rule all:      
    input:     
        expand("{outdir}/truvari/{sample}.{ref_name}.consensus_SVs_vs_GIAB/formatted_summary.txt",
                outdir=outdir,
                sample=samples,
               ref_name=["GRCh38", "GRCh37"]),
        expand("{outdir}/truvari/{sample}.{ref_name}.{tool}_vs_GIAB/formatted_summary.txt",
               outdir=outdir,
               sample=samples,
               ref_name=["GRCh37", "GRCh38"],
               tool=["sniffles", "cuteSV", "SVIM", "pbsv"]),
        expand("{outdir}/deep_variant/{sample}.{ref_name}.deep_variant.small_vars.vcf.gz",
               outdir=outdir,
               sample=samples,
               ref_name=["GRCh38", "GRCh37"]),
        expand("{outdir}/whatshap/{sample}.{ref_name}.whatshap.phasing_stats.txt",
                outdir=outdir,
                sample=samples,
                ref_name="GRCh38"),
        expand("{outdir}/alignment/{sample}.{ref_name}.sorted.whatshap_PHASED.bam",
                outdir=outdir,
                sample=samples,
                ref_name="GRCh38"),
        expand("{outdir}/whatshap/{sample}.{ref_name}.whatshap.compare_switch_errors.tsv",
               outdir=outdir,
               sample=samples,
               ref_name="GRCh38")
        expand("{outdir}/happy/{sample}.{ref_name}.deepvariant_vs_{truthset}/test.summary.csv",
               outdir=outdir,
               sample=samples,
               ref_name="GRCh38",
               truthset=["NISTv4.2", "CMGR"]),

rule align_minimap2:
    threads: 16
    resources:
        mem=128,
        time=32
    priority: 10
    input:
        fastq=lambda w: fastq_dict[w.sample],
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmpdir = lambda w: join(tmpdir, w.sample, w.ref_name),
    output:
        bam = join(outdir, "alignment/{sample}.{ref_name}.sorted.bam")
    conda: 'envs/LR_workflow.yaml'
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}

        cp {input.fastq} ./ 
        fastq=$(basename {input.fastq})

        cp {input.ref} ./
        ref=$(basename {input.ref})

        minimap2 -x map-hifi -a --MD -Y -t {threads} $ref $fastq | samtools view -hb > tmp.unsorted.bam
        samtools sort -m 1G tmp.unsorted.bam > {output.bam}
        rm -rf {params.tmpdir}
    """

rule downsample:
  threads: 1 
  resources:
      mem=24,
      time=4
  input:
    rules.align_minimap2.output[0]
  output:
    bam = join(outdir, "alignment/{sample}.{ref_name}.downsampled.bam"),
  conda: 'envs/LR_workflow'
  shell: """
    samtools view -s .01 -b {input} > {output.bam}
  """

rule downsample_nanoplot:
  threads: 1 
  resources:
    mem=24,
    time=6
  input:
      bam = rules.downsample.output[0],
      bai = rules.downsample.output[0] + '.bai'
  output:
      join(outdir, "nanoplot/downsampled-{sample}_{ref_name}/NanoPlot-report.html")
  conda: "envs/LR_workflow.yaml"
  shell: """
        out_dir=$(dirname {output})
        NanoPlot -t {threads} --raw --bam {input.bam} -o $out_dir
  """

rule index_bam:
    threads: 1
    resources:
        mem=32,
        time=8
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    conda: 'envs/LR_workflow.yaml'
    shell: """
        samtools index {input}
    """

rule nanoplot:
    threads: 4
    resources:
        mem=50,
        time=10
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai'
    output:
        join(outdir, "nanoplot/nanoplot-{sample}_{ref_name}/NanoPlot-report.html")
    conda: "envs/LR_workflow.yaml"
    shell: """
        out_dir=$(dirname {output})
        NanoPlot -t {threads} --bam {input.bam} --drop_outliers -o $out_dir
    """

rule cuteSV:
    threads: 4
    resources:
        mem=32,
        time=10
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"],
    params:
        tmpdir = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/cuteSV"),
        tmp_vcf = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/cuteSV/cuteSV_tmp.vcf")
    output:
        sv_vcf = join(outdir, "SV_calls/{sample}/{ref_name}/cuteSV/{sample}.{ref_name}.SVs.cuteSV.vcf")
    conda: "envs/LR_workflow.yaml"
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}
        ## Using the suggested advanced parameters on cuteSV for ONT data
        ## -s 3 --> decrease read support to call SVs to 3 given low-depth data
        ## -L 10000000 --> maximum SV length to call is 10Mb
        cuteSV -t {threads} \
            --sample {wildcards.sample}_cuteSV \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --report_readid \
            --min_support 5 \
            --min_size 30 \
            --max_size 10000000 \
            --genotype \
            {input.bam} {input.ref} {params.tmp_vcf} {params.tmpdir}

        # Sort VCF
        #bcftools sort {params.tmp_vcf} > {output.sv_vcf}
        mv {params.tmp_vcf} {output.sv_vcf}
        rm -rf {params.tmp_vcf}
    """

rule SVIM:
    threads: 1
    resources:
        mem=24,
        time=10
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"],
    params:
        tmpdir = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/SVIM"),
        tmp_vcf = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/SVIM/variants.vcf")
    output:
        sv_vcf = join(outdir, "SV_calls/{sample}/{ref_name}/SVIM/{sample}.{ref_name}.SVs.SVIM.vcf")
    conda: "envs/LR_workflow.yaml"
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}
        svim alignment \
            --sample {wildcards.sample}_SVIM \
            --minimum_depth 5 \
            --min_sv_size 30 \
            --max_sv_size 10000000 \
            --read_names \
            {params.tmpdir} {input.bam} {input.ref}

        ## SVIM reports all variants even low quality ones with little read support
        # Here we filter the full variant call set down to those with SUPPORT >=3
        awk '{{ if($1 ~ /^#/) {{ print $0 }} else {{ if($6>=3) {{ print $0 }} }} }}' {params.tmp_vcf} |
            sed 's/READS=/RNAMES=/g' | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' > {output.sv_vcf}
        rm -rf {params.tmp_vcf}
    """

rule pbsv:
    threads: 16
    resources:
      mem=240,
      time=24
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"],
    params:
        tmpdir = lambda w: join(tmpdir, w.sample, w.ref_name, "pbsv"),
        tmp_svsig = lambda w: join(tmpdir, w.sample, w.ref_name, "pbsv/pbsv_tmp.svsig.gz"),
        vcf = lambda w: join(tmpdir, w.sample, w.ref_name, "pbsv/pbsv_tmp.vcf"),
        trf_bed = lambda w: config["TRF_annotation"][w.ref_name]
    output:
        vcf = join(outdir, "SV_calls/{sample}/{ref_name}/pbsv/{sample}.{ref_name}.SVs.pbsv.vcf.gz"),
        tbi = join(outdir, "SV_calls/{sample}/{ref_name}/pbsv/{sample}.{ref_name}.SVs.pbsv.vcf.gz.tbi")
    conda: "envs/pbsv.yaml"
    shell: """
      ml bcftools
      mkdir -p {params.tmpdir} 
      cd {params.tmpdir}

      pbsv discover --hifi -s {wildcards.sample} --tandem-repeats {params.trf_bed} {input.bam} {params.tmp_svsig}
      pbsv call --hifi -j {threads} \
        --min-sv-length 30 \
        {input.ref} {params.tmp_svsig} {params.vcf}

      bcftools sort -Oz -o {output.vcf} {params.vcf}
      bcftools index --tbi {output.vcf}
      rm -rf {params.tmpdir}
    """


rule sniffles:
    threads: 4
    resources:
        mem=48,
        time=34
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"],
    params:
        tmpdir = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/sniffles"),
        tmp_vcf = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/sniffles/sniffles_tmp.vcf"),
        trf_bed = lambda w: config["TRF_annotation"][w.ref_name]
    output:
        vcf = join(outdir, "SV_calls/{sample}/{ref_name}/sniffles/{sample}.{ref_name}.SVs.sniffles.vcf"),
        snf = join(outdir, "SV_calls/{sample}/{ref_name}/sniffles/{sample}.{ref_name}.SVs.sniffles.snf")
    conda: "envs/LR_workflow.yaml"
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}
        sniffles \
            -t {threads} \
            --sample-id "{wildcards.sample}_sniffles" \
            --input {input.bam} \
            --minsvlen 30 \
            --minsupport 5 \
            --vcf {output.vcf} \
            --snf {output.snf} \
            --output-rnames \
            --tandem-repeats {params.trf_bed}
    """

rule jasmine_intrasample:
    threads: 1
    resources:
        mem=24,
        time=10
    priority: 10
    input:
        sniffles = rules.sniffles.output[0], 
        SVIM = rules.SVIM.output[0],
        cuteSV = rules.cuteSV.output[0], 
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmp_dir = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/jasmine"),
        vcf_list = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/jasmine/SV_vcfs.list"),
        tmp_vcf = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/jasmine/jasmine.tmp.vcf"),
        project_dir = config["project_dir"]
    output:
        temp(join(outdir, "SV_calls/{sample}/{ref_name}/{sample}.{ref_name}.consensus_SVs.jasmine.vcf")),
    conda: "envs/jasmine.yaml"
    shell: """
        mkdir -p {params.tmp_dir}
        echo {input.sniffles} > {params.vcf_list}
        echo {input.SVIM} >> {params.vcf_list}
        echo {input.cuteSV} >> {params.vcf_list}
        jasmine -Xmx18g --allow_intrasample --normalize_type --ignore_strand \
                        --file_list={params.vcf_list} --output_genotypes \
                        --dup_to_ins min_support=2 genome_file={input.ref} out_dir={params.tmp_dir} out_file={params.tmp_vcf}
        cat {params.tmp_vcf} | \
            python {params.project_dir}/scripts/consensus_genotype.py \
            > {output}
        rm -rf {params.tmp_dir}
    """ 
rule iris_refine:
    threads: 4
    resources:
        mem=32,
        time=24
    input:
        vcf = rules.jasmine_intrasample.output[0],
        bam = rules.align_minimap2.output[0],
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmp_dir = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/iris"),
        vcf_list = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/iris/SV_vcfs.list"),
        bam_list = lambda w: join(tmpdir, w.sample, w.ref_name, "SV_calls/iris/bams.list")
    output:
        temp(join(outdir, "SV_calls/{sample}/{ref_name}/{sample}.{ref_name}.consensus_SVs.jasmine_irisRefined.vcf"))
    conda: "envs/jasmine.yaml"
    shell: """
        mkdir -p {params.tmp_dir}
        cd {params.tmp_dir}
        ls {input.vcf} > {params.vcf_list}
        ls {input.bam} > {params.bam_list}
        jasmine -Xmx 24 --preprocess_only --run_iris iris_args=keep_long_variants,vcf_out={output} \
            genome_file={input.ref} file_list={params.vcf_list} bam_list={params.bam_list} \
            out_dir={params.tmp_dir} out_file={output} threads={threads}
        rm -rf {params.tmp_dir}
    """

rule sort_index_consensus_vcf:
  threads: 1 
  resources:
      mem=32,
      time=4
  input:
      vcf = rules.iris_refine.output[0]
  params:
      project_dir = config["project_dir"],
  output:
      vcf = join(outdir, "SV_calls/{sample}.{ref_name}.consensus_SVs.jasmine_irisRefined.sorted.vcf.gz"),
      tbi = join(outdir, "SV_calls/{sample}.{ref_name}.consensus_SVs.jasmine_irisRefined.sorted.vcf.gz.tbi")
  shell: """ 
    ml bcftools
    cat {input} | python {params.project_dir}/scripts/prep_for_sorting.py | bcftools sort -Oz -o {output.vcf}
    bcftools index --tbi {output.vcf}
  """

rule sort_index_vcf:
  threads: 1 
  resources:
      mem=32,
      time=4
  input:
      vcf = "{prefix}.vcf"
  output:
      vcf = "{prefix}.vcf.gz",
      tbi = "{prefix}.vcf.gz.tbi"
  shell: """ 
    ml bcftools
    # remove RNAMES from INFO field which causes issues for ONT duplex data
    cat {input} | sed 's/RNAMES=.*=.*;//g' | sed 's/;RNAMES=[^\t]*//g' | bcftools sort -Oz -o {output.vcf}

    bcftools index --tbi {output.vcf}
  """

rule truvari_single_tool:
  threads: 1 
  resources:
    mem=32,
    time=4
  input:
    vcf=join(outdir, "SV_calls/{sample}/{ref_name}/{tool}/{sample}.{ref_name}.SVs.{tool}.vcf.gz"),
    tbi=join(outdir, "SV_calls/{sample}/{ref_name}/{tool}/{sample}.{ref_name}.SVs.{tool}.vcf.gz.tbi")
  params:
    truthset=lambda w: config["SV_truth_set"][w.ref_name],
    truthbed = lambda w: config["SV_truth_bed"][w.ref_name]
  output:
    join(outdir, "truvari/{sample}.{ref_name}.{tool}_vs_GIAB/formatted_summary.txt")
  conda: "truvari"
  shell: """
    outdir=$(dirname {output})
    rm -rf $outdir
    truvari bench \
        -c  {input.vcf} \
        -b  {params.truthset} \
        --includebed {params.truthbed} \
        --sizemax 100000 --dup-to-ins -p 0 \
        -r 2000 --chunksize 2000 --no-ref a  \
        --sizemin 50 --sizefilt 50 \
        --passonly -o $outdir
    truthset=$(basename {params.truthset})
    cat $outdir/summary.json |  python3 -c "import sys,json; r=json.load(sys.stdin); print('{wildcards.sample}', '{wildcards.tool}', '$truthset', '\t'.join(list(map(lambda x: str(r[x]), ['TP-base','FP','FN','precision','recall','f1']))))" > {output}
  """

rule truvari_consensus:
  threads: 1 
  resources:
    mem=32,
    time=4
  input:
    vcf=rules.sort_index_consensus_vcf.output[0]
  params:
    truthset=lambda w: config["SV_truth_set"][w.ref_name],
    truthbed = lambda w: config["SV_truth_bed"][w.ref_name]
  output:
    join(outdir, "truvari/{sample}.{ref_name}.consensus_SVs_vs_GIAB/formatted_summary.txt")
  conda: "truvari"
  shell: """
    outdir=$(dirname {output})
    rm -rf $outdir
    truvari bench \
        -c  {input.vcf} \
        -b  {params.truthset} \
        --includebed {params.truthbed} \
        --sizemax 100000 --dup-to-ins -p 0 \
        -r 2000 --chunksize 2000 --no-ref a \
        --sizemin 50 --sizefilt 50 \
        --passonly -o $outdir
    truthset=$(basename {params.truthset})
    cat $outdir/summary.json |  python3 -c "import sys,json; r=json.load(sys.stdin); print('{wildcards.sample}', 'consensus_SVs', '$truthset', '\t'.join(list(map(lambda x: str(r[x]), ['TP-base','FP','FN','precision','recall','f1']))))" > {output}
  """

rule mosdepth:
    resources:
        mem=32,
        time=4,
    threads: 4
    input:
        bam= rules.align_minimap2.output[0],
        bai=rules.align_minimap2.output[0]+'.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        prefix_dir=join(outdir, "{sample}/mosdepth"),
        prefix="{sample}.{ref_name}",
        window_size = 1000
    output:
        join(outdir, "{sample}/mosdepth/{sample}.{ref_name}.mosdepth.region.dist.txt")
    conda: 'envs/cnvcalling.yaml'
    shell: '''
        mkdir -p {params.prefix_dir}
        prefix="{params.prefix_dir}/{params.prefix}"
        mosdepth --by {params.window_size} -n -t {threads} -f {input.ref} $prefix {input.bam}
    '''

rule spectre:
    resources:
        mem=24,
        time=6
    threads: 1
    input:
        mosdepth = rules.mosdepth.output[0],
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        sampleid = "{sample}.{ref_name}.spectre_CNV",
        project_dir = config["project_dir"],
        window_size = 1000
    output:
        join(outdir, "{sample}/SV_calls/spectre_CNV/{sample}.{ref_name}.spectre_CNV.vcf")
    conda: 'envs/cnvcalling.yaml'
    shell: """
        outdir=$(dirname {output})
        mosdepth_dir=$(dirname {input.mosdepth})
        mkdir -p $outdir

        cd $outdir
        python3 {params.project_dir}/scripts/Spectre/spectre.py CNVCaller \
            --bin-size {params.window_size} \
            --coverage $mosdepth_dir \
            --sample-id {params.sampleid} \
            --output-dir $outdir \
            --reference {input.ref} 
    """

rule vamos:
    resources:
        mem=72,
        time=8
    threads: 1
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai'
    params:
        vntr_motifs = config["vntr_motifs"],
        vamos=config["vamos_path"],
        project_dir = config["project_dir"],
        sample = '{sample}'
    output:
        vcf = join(outdir, "{sample}/SV_calls/vamos/{sample}.{ref_name}.vntr.vamos.vcf"),
        tsv = join(outdir, "{sample}/SV_calls/vamos/{sample}.{ref_name}.vntr.vamos.tsv")
    conda: 'envs/vamos.yaml'
    shell: """
        {params.vamos} --read \
            -t {threads} \
            -b {input.bam} \
            -r {params.vntr_motifs} \
            -s {params.sample} \
            -o {output.vcf}
        cat {output.vcf} | python3 {params.project_dir}/scripts/vntr_vcf_convert_tsv.py > {output.tsv}
    """

rule deepvariant:
    threads: 32
    resources:
        mem=84,
        time=36
    input:
        bam = rules.align_minimap2.output[0],
        bai = rules.align_minimap2.output[0] + '.bai',
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"],
    params:
        tmpdir=join(tmpdir,"{sample}_{ref_name}_deepvariant"),
        model=lambda w: model_map[w.sample],
        par_bed=lambda w: config["PAR_bed"][w.ref_name]
    output:
        vcf = join(outdir, "deep_variant/{sample}.{ref_name}.deep_variant.small_vars.vcf.gz"),
    singularity: "docker://google/deepvariant:1.6.0"
    shell: """
        mkdir -p {params.tmpdir}
        cp {input.bam} {params.tmpdir}/input.bam
        cp {input.bai} {params.tmpdir}/input.bam.bai 
        cp {input.ref} {params.tmpdir}/{wildcards.ref_name}.ref.fa 
        cp {input.ref}.fai {params.tmpdir}/{wildcards.ref_name}.ref.fa.fai 
        cp {params.par_bed} {params.tmpdir}/GRCh38_par.bed 
        /opt/deepvariant/bin/run_deepvariant \
          --model_type={params.model} \
          --ref={params.tmpdir}/{wildcards.ref_name}.ref.fa \
          --reads={params.tmpdir}/input.bam \
          --output_vcf={params.tmpdir}/output.vcf.gz \
          --num_shards={threads} \
          --haploid_contigs="chrX,chrY"  \
          --par_regions_bed={params.tmpdir}/GRCh38_par.bed \
          --dry_run=false 

        bcftools sort -Oz -o {output.vcf} {params.tmpdir}/output.vcf.gz
    """

rule happy:
  threads: 1 
  resources:
    mem=24,
    time=12
  input:
    vcf = rules.deepvariant.output[0],
    vcf = labmda w: config["smallvar_vcf"][w.ref_name][w.sample]
    ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
  params:
    truthset= lambda w: config["smallvar_truth_set"][w.ref_name][w.truthset],
    truthbed = lambda w: config["smallvar_truth_bed"][w.ref_name][w.truthset]
    extra_opt = lambda w: "-l chr20" if w.truthset == "NISTv4.2" else ""
  output:
    join(outdir, "happy/{sample}.{ref_name}.deepvariant_vs_{truthset}/test.summary.csv")
  conda: "happy"
  shell: """
    outdir=$(dirname {output})
    hap.py  \
        {params.truthset} \
        {input.vcf} \
        -r {input.ref} \
        -f {params.truthbed} \
        -o $outdir/test \
        {params.extra_opt}
  """

rule whatshap_phase:
    threads: 1
    resources:
        mem=42,
        time=16
    input:
        vcf = rules.deepvariant.output[0],
        bam = rules.align_minimap2.output[0],
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmpdir = join(tmpdir, "{sample}_{ref_name}_whatshap"),
        biallelic = join(tmpdir, "{sample}_{ref_name}_whatshap/tmp.{sample}.{ref_name}.biallelic.sites.vcf"),
        tmp_vcf = join(tmpdir, "{sample}_{ref_name}_whatshap/tmp.{sample}.{ref_name}.deep_variant.small_vars.whatshap_PHASED.vcf")
    output:
        phased_vcf = join(outdir, "whatshap/{sample}.{ref_name}.deep_variant.small_vars.whatshap_PHASED.vcf.gz"),
        tbi = join(outdir, "whatshap/{sample}.{ref_name}.deep_variant.small_vars.whatshap_PHASED.vcf.gz.tbi"),
    conda: 'envs/phasing.yaml'
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}

        bcftools view -m2 -M2 -f "PASS" {input.vcf} > {params.biallelic}
        whatshap phase --ignore-read-groups --indels \
            --reference {input.ref} \
            -o {params.tmp_vcf}  \
            {params.biallelic} {input.bam}

        # compress and index file
        bcftools sort -Oz -o {output.phased_vcf} {params.tmp_vcf}
        bcftools index --tbi {output.phased_vcf}

        rm -rf {params.tmpdir}
    """

rule whatshap_stats:
    threads: 1
    resources:
        mem=12,
        time=6
    input:
        rules.whatshap_phase.output[0]
    output:
        phased_stats = join(outdir, "whatshap/{sample}.{ref_name}.whatshap.phasing_stats.txt"),
        phased_bed = join(outdir, "whatshap/{sample}.{ref_name}.whatshap.phasing_blocks.bed"),
        gtf = join(outdir, "whatshap/{sample}.{ref_name}.whatshap.phasing_blocks.gtf")
    conda: 'envs/phasing.yaml'
    shell: """
        whatshap stats {input} --tsv={output.phased_stats} --block-list={output.phased_bed} --gtf={output.gtf}
    """

rule whatshap_compare:
    threads: 1 
    resources:
      mem=32,
      time=10
    input:
      rules.whatshap_phase.output[0]
    params:
      truth_vcf = lambda w: config["phasing_truth_set"][w.ref_name]
    output:
      join(outdir,"whatshap/{sample}.{ref_name}.whatshap.compare_switch_errors.tsv")
    conda: 'envs/phasing.yaml'
    shell: """
      whatshap compare --names truth,{wildcards.sample} --tsv-pairwise {output} --ignore-sample-name {params.truth_vcf} {input}
    """

rule haplotag_bam:
    threads: 1
    resources:
        mem=42,
        time=16
    input:
        bam = rules.align_minimap2.output[0],
        phased_vcf = rules.whatshap_phase.output[0],
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    output:
        phased_bam = join(outdir, "alignment/{sample}.{ref_name}.sorted.whatshap_PHASED.bam"),
        bai = join(outdir, "alignment/{sample}.{ref_name}.sorted.whatshap_PHASED.bam.bai")
    conda: 'envs/phasing.yaml'
    shell: """
        whatshap haplotag --reference {input.ref} \
            --ignore-read-groups \
            -o {output.phased_bam}  \
            {input.phased_vcf} {input.bam}
        samtools index {output.phased_bam}
    """
