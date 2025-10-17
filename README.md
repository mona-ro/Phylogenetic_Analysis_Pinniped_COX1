# Phylogenetic Analysis of Pinniped COX1 Genes

This project investigates the evolutionary relationships of pinnipeds, sirenians, and related mammals using COX1 gene sequences. The workflow integrates multiple sequence alignment (MSA), Neighbor-Joining (NJ) and Maximum Likelihood (ML) tree construction, nucleotide substitution model selection, and bootstrap analysis to infer phylogenetic patterns and divergence timing.

## Key Tasks and Skills
- Performed **multiple sequence alignment** of DNA sequences using `msa` (R).  
- Constructed **phylogenetic trees** with `phangorn` (NJ and ML methods).  
- Conducted **model testing and selection** for nucleotide substitution models (GTR+G+I). 
- Tested standard nucleotide substitution models (JC, F81, HKY, SYM, GTR), selected best-fit model using AICc.
- Executed **bootstrap analysis** (100 replicates) to evaluate branch support.  
- Generated **high-resolution tree visualizations** (PNG/PDF).  
- Interpreted evolutionary relationships, assessed **monophyly vs paraphyly**, and estimated divergence times.

## Tools
R, `msa`, `phangorn`, `ggplot2`, `Biostrings` (DNAStringSet)

## Repository Contents
- `Pinnipeds_initdata.txt` — initial COX1 sequences for pinniped species.  
- `Pinnipeds_newdata.txt` — additional sequences for updated analysis.  
- `phylogeny_analysis.R` — R script performing MSA, NJ/ML trees, model selection, and bootstrap analysis.  
- High-resolution tree images (`NJ_tree.png`, `ML_tree.png`) for quick visualization of results.  
- `new_msa_output.txt` — text file of the aligned sequences for reference.

## Highlights
- Strong bootstrap support confirms evolutionary relationships among pinnipeds, fur seals/sea lions, and walruses.  
- Independent marine adaptations observed: pinnipeds are closer to dogs, sirenians to elephants, indicating **convergent evolution**.  
- Reproducible workflow suitable for extending to other mammalian COX1 datasets.