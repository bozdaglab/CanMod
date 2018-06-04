# CanMod

CanMod is a computational pipeline to identify gene regulatory network in cancer. It consists of 6 main step

* Step 1: Identify regulators (i.e., miRNAs or transcription factors (TFs)) for each differentially expressed (DE) genes 
* Step 2: Cluster the regulators based on the similarity of their target genes
* Step 3: Cluster the target genes based on their GO-Term similarity
* Step 4: Generate modules using the regulators and the target gene clusters
* Step 5: Refine each modules so that the expression of regulators and target genes in each module is correlated
* Step 6: Merge modules which share a large amount of similar targets

Source code to run Canmod could be found in run_canmod.R

