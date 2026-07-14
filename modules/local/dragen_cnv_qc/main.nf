process DRAGEN_CNV_QC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-karyoploter:1.32.0--r44hdfd78af_0':
        'quay.io/biocontainers/bioconductor-karyoploter:1.32.0--r44hdfd78af_0' }"

    input:
    path(dragen_dirs)

    output:
    path("*.csv") , emit: metrics
    path("*.json"), emit: json
    path("*.png") , emit: plots
    path "versions.yml"
    tuple val("${task.process}"), val('karyoploter'), eval(''), topic: versions, emit: versions_karyoploter

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    analyze_cnv_metrics.R \\
        --dragen $dragen_dirs \\
        --prefix "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep "R version")
        karyoploteR: \$(Rscript -e "cat(as.character(utils::packageVersion('karyoploteR')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    """
    echo $args
    
    touch cnv-metrics.csv
    touch pur-cov-allmodel-stub.png
    touch pur-cov-maxmodel-stub.png
    touch karyoplot-stub-tn.png
    touch karyoplot-stub-tumor.baf.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | grep "R version")
        karyoploteR: \$(Rscript -e "cat(as.character(utils::packageVersion('karyoploteR')))")
    END_VERSIONS
    """
}
