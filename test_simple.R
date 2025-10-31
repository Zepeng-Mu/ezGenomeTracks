# Simple test script
library(ezGenomeTracks)
library(ggplot2)
library(dplyr)

# Test basic functionality
data(example_signal)
print("Data loaded successfully")

# Test parse_region
region_gr <- parse_region("chr1:1000000-2000000")
print("Region parsed successfully")
print(paste("Start:", GenomicRanges::start(region_gr)))
print(paste("End:", GenomicRanges::end(region_gr)))

# Test ez_signal
p <- ez_signal(example_signal, "chr1:1000000-2000000")
print("ez_signal created successfully")
print(class(p))
