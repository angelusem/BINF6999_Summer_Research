# Run from the directory containing this script and aligned_PCNA.fasta.
# install.packages('circlize')
library(circlize)
input <- 'circos_parent_only/aligned_PCNA_parent_only_sep30.fasta'
x <- readLines(input)
h <- which(startsWith(x, '>'))
ends <- c(h[-1]-1, length(x))
names <- sub('^>', '', x[h])
seqs <- vapply(seq_along(h), function(i) toupper(paste0(x[(h[i]+1):ends[i]], collapse='')), '')
stopifnot(length(unique(nchar(seqs))) == 1)
m <- do.call(rbind, strsplit(seqs, ''))
rownames(m) <- sub('_parent$', '', names)
# ordering totthed and baleens together
species_order <- c(
  # Toothed whales
  "BelugaWhale",
  "BottlenosDolphin",
  "Narwhal",
  "Orca",
  "PygmySpermWhale",
  "SpermWhale",
  # Baleen whales
  "BlueWhale",
  "NARightWhale",
  "Bowhead",
  "HumpbackWhale",
  "GrayWhale",
  "MinkeWhale"
)

stopifnot(setequal(species_order, rownames(m)))
m <- m[species_order, , drop = FALSE]
# proper names
display_names <- c(
  BelugaWhale     = "Beluga whale",
  BottlenosDolphin = "Bottlenose dolphin",
  Narwhal         = "Narwhal",
  Orca            = "Orca",
  PygmySpermWhale = "Pygmy sperm whale",
  SpermWhale      = "Sperm whale",
  BlueWhale       = "Blue whale",
  NARightWhale    = "North Atlantic right whale",
  Bowhead         = "Bowhead whale",
  HumpbackWhale   = "Humpback whale",
  GrayWhale       = "Gray whale",
  MinkeWhale      = "Minke whale"
)
group_colors <- c(
  "Toothed whales" = "#fc9272", # purple: #663399, dark red:#8B0000
  "Baleen whales"  = "#555555"
)

species_groups <- c(
  rep("Toothed whales", 6),
  rep("Baleen whales", 6)
)

label_colors <- unname(group_colors[species_groups])
rownames(m) <- unname(display_names[rownames(m)])
stopifnot(all(m %in% c('A','C','G','T'))) # This supplied alignment has no gaps/unknowns.
n <- nrow(m); L <- ncol(m)
variable <- apply(m, 2, function(z) length(unique(z)) > 1)
pos <- which(variable)
conservation <- apply(m, 2, function(z) max(table(z))/length(z))
cols <- c(A='#009E73', C='#0072B2', G='#E69F00', T='#CC79A7')
write.csv(data.frame(position=pos, t(m[,pos,drop=FALSE]), check.names=FALSE),
          'PCNA_variable_sites.csv', row.names=FALSE)
write.csv(data.frame(position=seq_len(L), conservation=conservation,
                     variable=variable), 'PCNA_conservation.csv', row.names=FALSE)

# One sector represents the shared CDS axis. Rings follow FASTA order,
# outermost to innermost. No clustering or reference species is used.
draw_circle <- function() {
  circos.clear()
  circos.par(start.degree=90, gap.degree=42, cell.padding=c(0,0,0,0),
             track.margin=c(0.002,0.002), points.overflow.warning=FALSE)
  circos.initialize('PCNA', xlim=c(0,L))
  circos.trackPlotRegion(ylim=c(0,1), track.height=0.10, bg.border=NA,
                         panel.fun=function(x,y) {
                           circos.rect(0,0,L,1,col='#F2F2F2',border=NA)
                           circos.lines(seq_len(L)-0.5,conservation,col='#333333',lwd=1)
                           circos.axis(h='top', major.at=c(0,seq(100,700,100),L),
                                       labels=c(1,seq(100,700,100),L), labels.cex=0.55)
                           circos.yaxis(side='left',at=c(0.5,1),labels=c('0.5','1'),labels.cex=0.5)
                         })
  for (i in seq_len(n)) {
    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = 0.043,
      bg.border = NA,
      panel.fun = function(x, y) {
        ring_background <- if (species_groups[i] == "Toothed whales") {
          "#fee0d2"  # Light purple #E6E6FA, light pink: "#FFB6C1"
        } else {
          "#E5E5E5"  # Light gray
        }
        
        circos.rect(0, 0, L, 1, col = ring_background, border = NA)
        circos.rect(
          pos - 1, 0, pos, 1,
          col = cols[m[i, pos]], border = NA
        )
        
        # Species labels colored by whale group
        theta <- get.cell.meta.data("cell.end.degree") * pi / 180
        radius <- mean(c(
          get.cell.meta.data("cell.top.radius"),
          get.cell.meta.data("cell.bottom.radius")
        ))
        
        text(
          radius * cos(theta),
          radius * sin(theta),
          rownames(m)[i],
          srt = get.cell.meta.data("cell.end.degree") - 90,
          adj = c(-0.08, 0.5),
          cex = 0.56,
          xpd = NA,
          col = label_colors[i]
        )
      }
    )
  }
  
  text(0, 0.10, "PCNA Coding Sequence", cex = 1.1, font = 2)
  text(
    0, -0.02,
    sprintf("%d species, %d bases\n%d variable sites", n, L, length(pos)),
    cex = 0.8
  )
  
  # Nucleotide legend
  legend(
    "bottomleft",
    legend = c(
      names(cols),
      "Toothed whales",
      "Baleen whales"
    ),
    fill = c(cols, "#fee0d2", "#E5E5E5"),
    border = NA,
    bty = "n",
    cex = 0.8
  )
  
  
  mtext(
    "Outer track represents the proportion of species carrying the most common nucleotide at a particular alignment position",
    side = 1, line = 0, cex = 0.7
  )
  
  circos.clear()
}
pdf('PCNA_circular.pdf',width=10,height=10); draw_circle(); dev.off()
png('PCNA_circular.png',width=2400,height=2400,res=240); draw_circle(); dev.off()

# Enlarged exact states at all variable columns: matching letters/colors
# in a column indicate shared states. Spacing here is intentionally uniform.
pdf('PCNA_variable_sites.pdf',width=12,height=6)
par(mar=c(6,12,4,2))
plot(NA,xlim=c(0.5,length(pos)+0.5),ylim=c(0.5,n+0.5),
     xlab='CDS nucleotide position (variable sites only)',ylab='',axes=FALSE,
     main='PCNA: shared nucleotide states at variable sites',xaxs='i',yaxs='i')
for(i in seq_len(n)) for(j in seq_along(pos)) {
  yy <- n-i+1
  rect(j-0.5,yy-0.5,j+0.5,yy+0.5,col=cols[m[i,pos[j]]],border='white')
  text(j,yy,m[i,pos[j]],cex=0.8)
}
axis(1,at=seq_along(pos),labels=pos,las=2)
axis(2,at=n:1,labels=rownames(m),las=2,tick=FALSE)
dev.off()
cat(sprintf('%d / %d columns conserved (%.2f%%).\n',sum(!variable),L,100*mean(!variable)))
