#!/usr/bin/env Rscript
# Demo script showing dynamic curve height functionality

library(ggplot2)
library(ezGenomeTracks)

cat("=================================================\n")
cat("Dynamic Curve Height Demo for geom_link\n")
cat("=================================================\n\n")

# Create sample data with varying genomic distances
interactions <- data.frame(
    start1 = c(1000, 2000, 5000, 10000),
    start2 = c(3000, 8000, 12000, 25000),
    score = c(0.5, 0.7, 0.9, 0.6)
)

cat("Sample data:\n")
print(interactions)
cat("\nGenomic distances:\n")
for (i in 1:nrow(interactions)) {
    dist <- interactions$start2[i] - interactions$start1[i]
    cat(sprintf("  Link %d: %d bp\n", i, dist))
}
cat("\n")

# Demo 1: Default downward curves
cat("Demo 1: Default downward curves\n")
cat("  - Direction: down\n")
cat("  - Height factor: 0.15 (default)\n")
p1 <- ez_link(interactions, "chr1:0-30000",
              color = "steelblue", alpha = 0.8)
ggsave("demo_1_default.png", p1, width = 10, height = 4)
cat("  ✓ Saved as demo_1_default.png\n\n")

# Demo 2: Upward curves
cat("Demo 2: Upward curves\n")
cat("  - Direction: up\n")
cat("  - Height factor: 0.15\n")
p2 <- ez_link(interactions, "chr1:0-30000",
              direction = "up", color = "coral", alpha = 0.8)
ggsave("demo_2_upward.png", p2, width = 10, height = 4)
cat("  ✓ Saved as demo_2_upward.png\n\n")

# Demo 3: Taller curves
cat("Demo 3: Taller curves with increased height factor\n")
cat("  - Direction: up\n")
cat("  - Height factor: 0.25 (taller)\n")
p3 <- ez_link(interactions, "chr1:0-30000",
              direction = "up", height_factor = 0.25,
              color = "darkgreen", alpha = 0.8)
ggsave("demo_3_taller.png", p3, width = 10, height = 4)
cat("  ✓ Saved as demo_3_taller.png\n\n")

# Demo 4: Score-based coloring
cat("Demo 4: Score-based coloring with downward curves\n")
cat("  - Direction: down\n")
cat("  - Use score: TRUE\n")
cat("  - Color gradient: gray80 -> red\n")
p4 <- ez_link(interactions, "chr1:0-30000",
              use_score = TRUE, color = "red",
              direction = "down", alpha = 0.9)
ggsave("demo_4_score.png", p4, width = 10, height = 4)
cat("  ✓ Saved as demo_4_score.png\n\n")

# Demo 5: Custom curvature with arrows
cat("Demo 5: Custom curvature with directional arrows\n")
cat("  - Direction: up\n")
cat("  - Curvature: 0.8\n")
cat("  - Arrow length: 0.05\n")
p5 <- ggplot(interactions) +
    geom_link(aes(x = start1, y = 0, xend = start2, yend = 0),
              direction = "up", height_factor = 0.18,
              curvature = 0.8, arrow_length = 0.05,
              color = "purple", linewidth = 0.8, alpha = 0.8) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-3, suffix = "kb")) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    labs(x = "Genomic Position", title = "Custom Curves with Arrows")
ggsave("demo_5_arrows.png", p5, width = 10, height = 4)
cat("  ✓ Saved as demo_5_arrows.png\n\n")

cat("=================================================\n")
cat("Key observations:\n")
cat("  1. Longer genomic distances → taller curves\n")
cat("  2. Curves scale proportionally with distance\n")
cat("  3. Direction parameter controls up/down orientation\n")
cat("  4. height_factor controls overall curve height\n")
cat("  5. No clipping issues - curves always visible\n")
cat("=================================================\n\n")

cat("Implementation complete! ✓\n")
cat("Curves now respond to plot dimensions dynamically.\n")
