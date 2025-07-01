---
title: 'Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assays'
tags:
  - R
  - bioinformatics
  - antigenic cartography
  - virology
authors:
  - name: Omid Arhami
    orcid: 0009-0005-2681-6598
    affiliation: 1, 3
  - name: Pejman Rohani
    orcid: 0000-0002-7221-3801
    affiliation: 1, 2, 3
affiliations:
 - name: Institute of Bioinformatics, University of Georgia, Athens, GA 30602, USA
   index: 1
 - name: Odum School of Ecology, University of Georgia, Athens, GA 30602, USA
   index: 2
 - name: Center for Influenza Disease \& Emergence Research (CIDER), University of Georgia, Athens, GA 30602, USA
   index: 3
date: 2025
bibliography: paper.bib
---

@article{10.1093/bioinformatics/btaf372,
    author = {Arhami, Omid and Rohani, Pejman},
    title = {Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assays},
    journal = {Bioinformatics},
    pages = {btaf372},
    year = {2025},
    month = {06},
    abstract = {Understanding antigenic evolution through cross-reactivity assays is crucial for tracking rapidly evolving pathogens requiring regular vaccine updates. However, existing cartography methods, commonly based on multidimensional scaling (MDS), face significant challenges with sparse and complex data, producing incomplete and inconsistent maps. There is an urgent need for robust computational methods that can accurately map antigenic relationships from incomplete experimental data while maintaining biological relevance, especially given that more than 95\% of possible measurements could be missing in large-scale studies.We present Topolow, an algorithm that transforms cross-reactivity and binding affinity measurements into accurate positions in a phenotype space. Using a physics-inspired model, Topolow achieved comparable prediction accuracy to MDS for H3N2 influenza and 56\% and 41\% improved accuracy for Dengue and HIV, while maintaining complete positioning of all antigens. The method effectively reduces experimental noise and bias, determines optimal dimensionality through likelihood-based estimation, avoiding distortions due to insufficient dimensions, and demonstrates orders of magnitude better stability across multiple runs. We also introduce antigenic velocity vectors, which measure the rate of antigenic advancement of each isolate per unit of time against its temporal and evolutionary related background, revealing the underlying antigenic relationships and cluster transitions.Topolow is implemented in R and freely available at [https://doi.org/10.5281/zenodo.15620983] and [https://github.com/omid-arhami/topolow].Available at Bioinformatics online.},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btaf372},
    url = {https://doi.org/10.1093/bioinformatics/btaf372},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btaf372/63582086/btaf372.pdf},
}






# Summary

Antigenic cartography is a crucial tool for visualizing and analyzing antigenic relationships between viruses. Here we present `topolow`, an R package implementing a novel physics-inspired algorithm for antigenic mapping. The algorithm uses spring forces and repulsive interactions to optimize point configurations in high-dimensional spaces while handling missing and thresholded measurements.
