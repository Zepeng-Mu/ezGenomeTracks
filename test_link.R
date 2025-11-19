library(ggplot2)
library(ezGenomeTracks)

# Create dummy data
df <- data.frame(
    start1 = c(1000, 2000, 3000),
    end1 = c(1100, 2100, 3100),
    start2 = c(5000, 6000, 7000),
    end2 = c(5100, 6100, 7100),
    score = c(0.8, 0.6, 0.9)
)

# Test geom_link directly
p1 <- ggplot(df) +
    geom_link(aes(x = start1, y = 0, xend = start2, yend = 0))
print(p1)

# Test ez_link
p2 <- ez_link(df, "chr1:1-10000", color = "blue")
print(p2)

# Test ez_link with score
p3 <- ez_link(df, "chr1:1-10000", use_score = TRUE, color = "red")
print(p3)

print("All tests passed!")
