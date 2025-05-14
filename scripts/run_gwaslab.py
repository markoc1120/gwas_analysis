import gwaslab as gl

mysumstats_linear = gl.Sumstats(
    "gwas/linear_clean.tsv",
    snpid="SNP",
    chrom="CHR",
    pos="BP",
    ea="A1",
    beta="BETA",
    p="P",
    n="NMISS",
    build="19",
)

mysumstats_linear.plot_mqq(mode="qqm", anno=True, save="plots/linear_gwaslab_mqq.png")
mysumstats_linear.plot_mqq(
    mode="b", bwindowsizekb=500, save="plots/linear_gwaslab_manhattan.png"
)

mysumstats_mlma = gl.Sumstats(
    "gwas/mlma_clean.tsv",
    snpid="SNP",
    chrom="CHR",
    pos="BP",
    ea="A1",
    nea="A2",
    eaf="FREQ",
    beta="B",
    se="SE",
    p="P",
    build="19",
)

mysumstats_mlma.plot_mqq(mode="qqm", anno=True, save="plots/mlma_gwaslab_mqq.png")
mysumstats_mlma.plot_mqq(
    mode="b", bwindowsizekb=500, save="plots/mlma_gwaslab_manhattan.png"
)
