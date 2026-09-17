#!/usr/bin/env Rscript
# plot_ggtree.R — publication figure from a tree bundle.
#
# Consumes a bundle directory produced by make_tree_bundle.py: tree.newick plus the tidy
# ann.annotations.tsv whose ID match rate the bundle manifest has already verified. Nothing here
# re-derives annotations; if the join is incomplete the script says so and stops, because a
# silently half-joined tree is the failure this whole bundle arrangement exists to prevent.
#
# Usage: Rscript plot_ggtree.R <bundle-dir> <out.pdf> [--clade TIP] [--circular]

suppressPackageStartupMessages({library(ggtree); library(treeio); library(ggplot2); library(dplyr)})
args <- commandArgs(trailingOnly=TRUE)
bundle <- args[1]; out <- args[2]
clade  <- if ("--clade" %in% args) args[which(args=="--clade")+1] else NA
circ   <- "--circular" %in% args

tr  <- read.newick(file.path(bundle,"tree.newick"), node.label="support")
ann <- read.delim(file.path(bundle,"ann.annotations.tsv"), stringsAsFactors=FALSE, na.strings="NA")

tips <- tr@phylo$tip.label
matched <- sum(tips %in% ann$label)
cat(sprintf("join: %d/%d tips matched the annotation table\n", matched, length(tips)))
if (matched < length(tips)) stop("annotation table does not cover every tip - wrong bundle?")

# optional: restrict to the smallest clade containing a named tip plus N others
if (!is.na(clade)) {
  nd <- which(tips == clade)
  if (!length(nd)) stop("tip not found: ", clade)
  anc <- ancestor(tr@phylo, nd)
  pick <- NA
  for (a in rev(anc)) { n <- length(offspring(tr@phylo, a, tiponly=TRUE)); if (n >= 20) { pick <- a; break } }
  if (!is.na(pick)) { tr <- treeio::tree_subset(tr, node=pick, levels_back=0); tips <- tr@phylo$tip.label
                      cat(sprintf("subset to %d tips around %s\n", length(tips), clade)) }
}

LIGCOL <- c("OA"="#C15C3E","TA"="#E0A458","TA/OA"="#D4813F","NE"="#8E3E29","DA"="#1A2C7A",
            "5-HT"="#3C8DAD","5-HT?"="#9EC5D4","ACh"="#7A7A7A","MEL"="#6B4E8E","HA"="#4F8A5B")
ann <- ann %>% mutate(
  lab  = ifelse(kind=="scan", paste0(gene,"  ", species), paste0(gene,"  ", species)),
  d332 = ifelse(is.na(res332), "", res332),
  bad  = !is.na(res332) & res332 != "D",
  subj = label %in% c("cel_dop-6","cel_dop-5") | grepl("^TEST_", label))

p <- ggtree(tr, size=.32, layout=if (circ) "circular" else "rectangular") %<+% ann +
  # support only where it is meaningful
  geom_text2(aes(subset = !is.na(as.numeric(sub("/.*","",support))) &
                          as.numeric(sub(".*/","",support)) >= 70,
                 label = sub(".*/","",support)),
             hjust=1.25, vjust=-.45, size=1.5, colour="grey45") +
  geom_tippoint(aes(colour=ligand), size=1.25, na.rm=TRUE) +
  scale_colour_manual(values=LIGCOL, na.value="transparent", name="ligand (measured)") +
  geom_tiplab(aes(label=lab, fontface=ifelse(subj,"bold","plain")), size=1.55, offset=.045) +
  geom_tiplab(aes(label=d332, colour=NULL), size=1.5, offset=.012,
              fontface="bold", family="mono", show.legend=FALSE) +
  geom_tiplab(aes(label=ifelse(is.na(coupling),"",coupling)), size=1.4,
              offset=.9, colour="grey40", show.legend=FALSE) +
  theme_tree2(legend.position=if (circ) "right" else c(.08,.88)) +
  theme(legend.key.size=unit(3.2,"mm"), legend.text=element_text(size=6),
        legend.title=element_text(size=6.5))
p <- p + xlim(0, max(p$data$x, na.rm=TRUE)*1.75)

h <- if (circ) 11 else max(4, length(tips)*0.075)
ggsave(out, p, width=if (circ) 11 else 9.5, height=h, limitsize=FALSE, device=cairo_pdf)
ggsave(sub("\\.pdf$",".png",out), p, width=if (circ) 11 else 9.5, height=h, dpi=320, limitsize=FALSE)
cat("->", out, "\n")
