#library("ape")
#library("Biostrings")
library("ggplot2")
library("ggtree")
library("optparse")

###### parser
option_list = list(
    make_option( c("-t","--file_tree"), type="character", default=NULL,
                 help="Input Newick tree file", metavar="character"),
    make_option( c("-a","--file_annotate"), type="character", default=NULL,
                 help="Input annotation .csv two colums list, separated by \'\\t\', labels should be:
                        1: Annotated only
                        2: Find only
                        3: Find & Annotated", metavar="character"),
    make_option( c("-s","--font_size"), type="double", default=1,
                 help="label font size, default=2", metavar="double"),
    make_option( c("-o","--file_output"), type="character", default="Rplots.pdf",
                 help="Output file name, default=\"Rplots.pdf\", default=2", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);



###### read tree and annotations
tree   <- read.tree(opt$file_tree)
tipcat <- read.csv(opt$file_annotate, sep="\t", col.name = c("allele","cat"), header = FALSE, stringsAsFactors = FALSE)

#plot tree
p <- ggtree(tree)#, layout="circular")
#plot labels
p <- p %<+% tipcat + geom_tiplab(aes(fill = factor(cat)), color = "black", size=opt$font_size, geom = "label", label.padding = unit(0.05,"lines"), label.size=0)
#adjust label colors
p <- p + scale_fill_manual(values=c('#FF7E7E','#A6E8FF','#D8A6FF'), 
                           label = c("Annotated only", "Find only","Find & Annotated"))
#adjust plot size to make the labels fit in
p <- p + xlim(0,1)
#adjust the legends
p <- p + theme(legend.position = c(0.8, 0.15),
               legend.title = element_blank(), legend.key= element_blank())
p <- p + geom_treescale(fontsize=4)


###### print pdf
pdf(opt$file_output)
print(p)

