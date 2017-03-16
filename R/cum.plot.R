#' @title Cumulative frequency plot
#' @author Peter Utnes \email{utnesp@gmail.com}
#' @param file.or.object Provide a .txt file or object where row.names are ensembl_gene_identifiers
#' @param legend.position Position of legend within plot. Defaults to c(0.8, 0.7) 
#' @param incl.mito.genes Include mitochondrial genes in plot. Defaults to FALSE.
#' @param cum.freq.out.file.txt If specified, writes out the cumulative frequencies to a file.
#' @param file.pdf  Where to store cum.plot. Defaults to getwd() and to filename "CPM.cumfreq.pdf".
#' @param file.html Where to store cum.plot. Defaults to getwd() and to filename "CPM.cumfreq.html".
#' @param color.numbers.to.use Specify colors to use for individual biotypes from the color palette brewer.pal.info[brewer.pal.info$category == 'qual',]. 
#' @param biotype.cutoff How many biotypes to represent in the plot. Defaults to 10.
#' @export
#' @import biomaRt
#' @import RColorBrewer
#' @import easybiomart
#' @import ggplot2
#' @import ggthemes
#' @import plotly
#' @import htmlwidgets
#' @examples See https://github.com/utnesp/cumulative.plot.R

cum.plot <- function(file.or.object, legend.position = c(0.8, 0.7), incl.mito.genes = F, cum.freq.out.file.txt = "", file.pdf = "", file.html = "", color.numbers.to.use = "", biotype.cutoff = "") {

    if (grepl(".txt", file.or.object) == TRUE) t <- read.delim(file.or.object, header = T)
    if ((class(file.or.object) == "matrix") == TRUE) t <- data.frame(file.or.object)
    t$CPM <- rowMeans(t)
    t$ensembl_gene_id <- row.names(t)
    t <- t[grep("ENSG", t$ensembl_gene_id), ]
    x <- ensg2ext_name_biotype(t$ensembl_gene_id)
    t <- merge(x, t, by = "ensembl_gene_id")
    x <- ensg2chr_start_end(t$ensembl_gene_id)
    t <- merge(x, t, by = "ensembl_gene_id")
    t$length <- abs(t$end_position) - abs(t$start_position)
    m <- table(t$gene_biotype)
    m <- data.frame(m)
    m <- m[order(-m$Freq), ]
    print(m)
    
    if (biotype.cutoff == "") {
        print("You have not set cutoff for which biotypes to include")
        print("Automatically using max ten biotypes")
        if (nrow(m) > 10) m <- m[1:10, ]
    } else {
        m <- m[m$Freq > biotype.cutoff, ]
    }
    
    print(m)
    t$gene_biotype <- ifelse(t$gene_biotype %in% m$Var1, t$gene_biotype, "other")
    
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    if (color.numbers.to.use == "") {
        print("Selecting colors automatically")
        colors <- col_vector[c(1,2,3,5,6,8,10, 14, 15, 17, 21, 22:60)]
        colors <- colors[1:(nrow(m)+1)] 
        pie(rep(1,nrow(m)+1), col=colors)
    } else {
        par(mfrow = c(2,1))
        print("You may choose from these colors:")
        pie(rep(1,length(col_vector)), col=col_vector)
        print("Selecting colors manually")
        colors <- col_vector[colors] 
        pie(rep(1,nrow(m)), col=colors)
        par(mfrow = c(1,1))
    }
    
    t$CPM_ensg <- paste(t$CPM, t$ensembl_gene_id, t$length, sep = "_")
    t <- melt(t, "gene_biotype", "CPM_ensg")
    t <- cSplit(t, "value", "_") 
    
    m <- table(t$gene_biotype)
    m <- data.frame(m)
    m <- m[order(-m$Freq), ]
    m$col.nr <- row.names(m)
    m$col <- factor(colors)
    colnames(m) <- c("gene_biotype", "bio_tot", "col.nr", "col")
    
    t <- merge(t, m, by = "gene_biotype")
    t$col.nr <- NULL; t$variable <- NULL
    colnames(t) <- c('gene_biotype','CPM','ensembl_gene_id','GeneLength','bio_tot','col')
    
    k <- aggregate(CPM ~ gene_biotype, t, sum)
    t <- merge(t, k, by = "gene_biotype", all = T)
    colnames(t) <- c('gene_biotype','CPM','ensembl_gene_id','GeneLength','bio_tot','col','sum')
    t$freq <- t$CPM / t$sum
    
    x <- ensg2ext_name(t$ensembl_gene_id)
    t <- merge(x, t, by = "ensembl_gene_id")
    t <- with(t, t[order(gene_biotype, freq), ])
    f <- setDT(t)
    e <- f[, CumulativeFrequency := cumsum(freq), by=list(gene_biotype)]
    e <- data.frame(e)
    e$gene_biotype <- paste(e$gene_biotype, " (", e$bio_tot, ")", sep = "")
    
    
    if (incl.mito.genes == T) {
        e.sub <- e
    } else {
        e.sub <- e[grep("MT-", e$external_gene_name, invert = T), ]
    }
    
    shapes = c(15, 16, 17, 0, 18, 3, 5, 2, 8, 24, 1, 25, 4, 7, 8, 9, 10, 26:30)
    
    s <- data.frame(gene_biotype = e.sub$gene_biotype[!duplicated(e.sub$gene_biotype)], shape = shapes[1:length(colors)] , color = colors)
    # e.sub$shape <- NULL; e.sub$color <- NULL
    e.sub <- merge(e.sub, s, by = "gene_biotype", all = T)
    
    # head(e.sub[grep("protein_coding", e.sub$gene_biotype), ])
    
    if (cum.freq.out.file.txt != "") {
        write.delim(e.sub, cum.freq.out.file.txt)
    }
    
    if (file.pdf == "") file.pdf <- paste(getwd(), "/CPM.cumfreq.pdf", sep = "")
    
    # p1 <- ggplot(e.sub, aes(CPM, CumulativeFrequency)) + scale_shape_identity(name = "Biotype") + geom_point(aes(color = factor(gene_biotype), shape = shape, text = paste("GENE: ", paste(external_gene_name, "(", ensembl_gene_id, ")", sep = ""), "LENGTH:", GeneLength)), size = 0.5) + theme_classic() + labs(x = "log2(CPM)", y = "Cumulative Frequency")
    # p1 <- ggplot(e.sub, aes(CPM, CumulativeFrequency)) + geom_point(aes(color = factor(gene_biotype), shape = factor(shape), text = paste("GENE: ", paste(external_gene_name, "(", ensembl_gene_id, ")", sep = ""), "LENGTH:", GeneLength)), size = 1) + theme_classic() + theme(legend.title = element_text(colour = "white")) + labs(x = "CPM") + guides(shape=FALSE) + labs(x = "log2(CPM)", y = "Cumulative Frequency")
    
    p1 <- ggplot(e.sub, aes(CPM, CumulativeFrequency)) + geom_point(aes(shape = gene_biotype, color = gene_biotype, text = paste("GENE: ", paste(external_gene_name, "(", ensembl_gene_id, ")", sep = ""), "LENGTH:", GeneLength)), size = 1) + theme_classic() + labs(x = "log2(CPM)", y = "Cumulative Frequency", color='Biotype') +
        scale_colour_manual(name = "Biotype",
                            labels = s$gene_biotype,
                            values = as.character(s$color))  +
        scale_shape_manual(name = "Biotype",
                           labels = s$gene_biotype,
                           values = s$shape) + 
        if(legend.position != "") theme(legend.text=element_text(size=10), legend.position = legend.position) + 
        if(legend.position == "") theme(legend.text=element_text(size=10), legend.justification = c(1,1), legend.position = c(1,1)) + 
        guides(shape=guide_legend(override.aes=list(size=5)))
    plot.new()
    p1
    ggsave(file.pdf, device = "pdf", dpi = 300)
    
    p2 <- ggplot(e.sub, aes(x = factor(1), fill = factor(gene_biotype))) + geom_bar(width = 1, show.legend = FALSE) + 
        scale_fill_manual(name = "Biotype",
                          labels = s$gene_biotype,
                          values = as.character(s$color))
    p2 <- p2 + coord_polar(theta = "y") + theme_classic() + labs(x="", y = "") + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) 
    p2
    dev.off()
    
    openPDF(file.pdf)
    system(paste("open", dirname(file.pdf)))
    
    p1 <- ggplotly(p1)
    if (file.html == "") file.html <- paste(getwd(), "/CPM.cumfreq.html", sep = "")
    htmlwidgets::saveWidget(p1, file.html, selfcontained = T, libdir = NULL)
    system(paste("open", file.html))
}