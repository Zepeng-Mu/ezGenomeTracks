---
name: user-test
description: Test the package with example data from human created examples and plotgardener data
author: "Zepeng Mu"
tags: [R, plotting, genomics, visualization]
---

This skill tests the core functionality of the package ezGenomeTracks using plotgardener data. The goal is to use real-world data instead of unit tests to reflect more realistic use cases for the package.

This test consists of two parts:

## Part 1: Use example data from package `plotgardenerData`. The example usage can be found here: https://phanstiellab.github.io/plotgardener/.

Test each of the following track types with example data:
- Gene track
- Feature track
- Coverage track
- Hi-C track
- Loop track

Also test the integration of multiple tracks in a single plot to ensure that they can be combined effectively without errors.

## Part 2: Use human created examples in the vignette of ezGenomeTracks here: "../../../vignettes/articles/".

In this part, make sure all vignettes can be compiled without errors and that the resulting plots are generated correctly. This includes testing any custom track types or configurations that are demonstrated in the vignettes.


Each time this skill is used, create an Rmd or quartet markdown file that includes code and output plot to generate each of the above track types using ezGenomeTracks. Ensure that the code runs without errors and produces the expected visualizations. Generate an HTML file from the Rmd or quartet markdown file to let human visually inspect the results.
