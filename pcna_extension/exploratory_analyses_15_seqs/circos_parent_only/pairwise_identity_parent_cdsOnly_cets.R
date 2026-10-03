library(ggplot2)

species <- rownames(identity)

# Extract each unique pair once
idx <- which(
  lower.tri(identity) & is.finite(identity),
  arr.ind = TRUE
)

heatmap_data <- data.frame(
  row_species = species[idx[, 1]],
  col_species = species[idx[, 2]],
  identity_pct = identity[idx]
)

# Preserve the same species order on both axes
heatmap_data$col_species <- factor(
  heatmap_data$col_species, levels = species
)
heatmap_data$row_species <- factor(
  heatmap_data$row_species, levels = rev(species)
)

clean_labels <- function(x) {
  x <- sub("_parent$", "", x)
  x <- gsub("([a-z])([A-Z])", "\\1 \\2", x)
  x <- sub("^NARight", "North Atlantic right", x)
  x[x == "Bottlenos Dolphin"] <- "Bottlenose Dolphin"
  x[x == "Bowhead"] <- "Bowhead Whale"
  x
}

p <- ggplot(
  heatmap_data,
  aes(col_species, row_species, fill = identity_pct)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", identity_pct)),
    color = "white",
    size = 5
  ) +
  scale_fill_gradient(
    low = "#40E0D0",     
    high = "#0000CD",    
    limits = c(98.70, 100),
    name = "Identity (%)"
  ) +
  scale_x_discrete(
    limits = species[-length(species)],
    labels = clean_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(species[-1]),
    labels = clean_labels,
    drop = FALSE
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 16,
      color = "black",
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(
      hjust = 1,
      size = 16,
      color = "black",
      margin = margin(r = 8)
    ),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.position = "right",
    plot.margin = margin(15, 20, 15, 15)
  )
  
p

ggsave(
  "pairwiseIdentity_cetCDS_heatmap.png",
  plot = p,
  width = 14,
  height = 13,
  units = "in",
  bg = "white"
)
  # Convert the color name "steelblue" to its hex code
  rgb(t(col2rgb("lavender")), maxColorValue = 255)
  # Output: "#4682B4"
  # light blue:"#BDD7E7, darker: '#0072B2 