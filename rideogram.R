require(RIdeogram)
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")

ideogram(karyotype = human_karyotype)


ideogram(karyotype, overlaid = NULL, label = NULL, label_type = NULL, synteny = NULL, colorset1, colorset2, width, Lx, Ly, output = "chromosome.svg")
ideogram(karyotype = human_karyotype, overlaid = gene_density)

convertSVG("chromosome.svg", device = "png")
