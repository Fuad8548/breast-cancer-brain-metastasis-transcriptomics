




# Repository structure

```
singlecell-brain-metastasis/
├── README.md                                      # Project documentation & key findings
├── config.R                                       # Global directory & 8GB RAM parameters
├── data/
│   ├── raw/                                       # GSE268662 10x triplets (git-ignored)
│   ├── signatures/                                # GSE25066_limma_DEGs.csv (2,938 probes)
│   └── serialized/                                # .rds checkpoints (01, 03, 04, 05)
├── scripts/
│   ├── 01_data_ingestion.R                        # Per-sample QC & 8GB RAM balancing
│   ├── 02_qc_report.R                             # Violin plots & per-sample audit table
│   ├── 03_normalization_clustering.R              # LogNormalize, Harmony & 17-cluster UMAP
│   ├── 04_chemoresistance_projection.R            # hgu133a.db probe bridge & AddModuleScore
│   ├── 04b_annotate_and_fix_dotplot.R             # Biological cell-type labels & clean DotPlot
│   └── 05_cellchat_neuro_immune.R                 # CellChat ligand-receptor crosstalk
└── outputs/
    ├── figures/
    │   ├── 02_qc_postfilter_violin.png
    │   ├── 03_umap_clusters_and_tissue.png
    │   ├── 04_lineage_markers_dotplot_clean.png
    │   ├── 04_annotated_celltypes_umap.png
    │   ├── 04_chemoresistance_projection_umap.png
    │   ├── 04_top_ml_drivers_featureplot.png
    │   ├── 05_cellchat_chemoresistant_circle.png
    │   └── 05_cellchat_ligand_receptor_bubble.png
    └── tables/
        ├── 02_qc_sample_audit.csv
        ├── 04_GSE25066_translated_gene_bridge.csv
        ├── 04_cluster_chemoresistance_scores.csv
        ├── 05_cellchat_all_interactions.csv
        └── 05_chemoresistant_tumor_to_immune_signaling.csv
```











































