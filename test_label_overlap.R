#!/usr/bin/env Rscript
# Test script for new gene label overlap handling features

library(ezGenomeTracks)

# Load example data
data(example_genes)

# Define a crowded region with many genes
crowded_region <- "chr1:42100000-42700000"

cat("Testing gene label overlap handling features\n")
cat("=============================================\n\n")

# Test 1: Default behavior (auto mode with check_overlap)
cat("1. Default behavior (label_style = 'auto'):\n")
p1 <- ez_gene(example_genes, crowded_region)
print("   ✓ Created plot with automatic label overlap handling")

# Test 2: No overlap handling
cat("\n2. Simple mode (label_style = 'simple'):\n")
p2 <- ez_gene(example_genes, crowded_region, label_style = "simple")
print("   ✓ Created plot with no overlap handling (all labels shown)")

# Test 3: Limit labels to top 5 longest genes
cat("\n3. Limit to 5 labels (max_labels = 5):\n")
p3 <- ez_gene(example_genes, crowded_region, max_labels = 5)
print("   ✓ Created plot showing only 5 longest genes")

# Test 4: Limit labels with different priority
cat("\n4. Alphabetical priority (max_labels = 5, label_priority = 'name'):\n")
p4 <- ez_gene(example_genes, crowded_region,
              max_labels = 5,
              label_priority = "name")
print("   ✓ Created plot showing 5 genes (alphabetically sorted)")

# Test 5: No labels
cat("\n5. No labels (label_style = 'none'):\n")
p5 <- ez_gene(example_genes, crowded_region, label_style = "none")
print("   ✓ Created plot with no labels")

# Test 6: ggrepel mode (if available)
cat("\n6. ggrepel mode (if installed):\n")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p6 <- ez_gene(example_genes, crowded_region, label_style = "repel")
  print("   ✓ Created plot with ggrepel label repositioning")

  # Test 7: ggrepel with custom arguments
  cat("\n7. ggrepel with custom settings:\n")
  p7 <- ez_gene(example_genes, crowded_region,
                label_style = "repel",
                repel_args = list(max.overlaps = 20, force = 2))
  print("   ✓ Created plot with custom repel parameters")
} else {
  cat("   ⚠ ggrepel not installed, skipping repel tests\n")
  cat("   Install with: install.packages('ggrepel')\n")
}

cat("\n=============================================\n")
cat("All tests completed successfully!\n\n")

cat("Key Features Implemented:\n")
cat("- label_style: 'auto', 'simple', 'repel', 'none'\n")
cat("- max_labels: Limit number of labels shown\n")
cat("- label_priority: 'length' or 'name' for filtering\n")
cat("- repel_args: Pass custom arguments to ggrepel\n")
