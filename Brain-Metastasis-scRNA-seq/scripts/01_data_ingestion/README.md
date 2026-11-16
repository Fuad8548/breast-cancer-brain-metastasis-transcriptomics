# Step-by-Step Technical Audit of Milestone 1 (01_data_ingestion.R)

1. Automated Archive Extraction & Exact 10x Triplet Pairing (Steps 1 & 2):
	- Unpacks GSE268662_RAW.tar automatically if .mtx.gz files aren't already extracted.
	- Parses exact sample prefixes (e.g., GSM8296168_DN_) to pair each sparse count matrix (matrix.mtx.gz) with its exact barcodes.tsv.gz and features.tsv.gz without fragile regex matching.
	- Initializes individual Seurat objects (min.cells = MIN_CELLS, min.features = 100) and tags every cell with its patient origin (`obj$sample_id <- sample_id`).
2. Per-Sample Quality Control & Balanced Downsampling (Step 3):
	- Dead/Dying Cell & Doublet Filtering: Computes mitochondrial transcript percentage (percent.mt via ^MT-) and removes broken cells (`nFeature_RNA <= MIN_FEATURES`), dying/lysed cells (`percent.mt >= MAX_MT_PCT`), and putative cell doublets (`nFeature_RNA >= MAX_FEATURES`).
	- Stratified Sample Balancing (MAX_CELLS_PER_SAMPLE, set.seed(42)): Caps the maximum number of cells per patient tumor sample.
		- Computational benefit: Keeps memory footprint strictly within 8GB RAM (`gc()` explicitly frees RAM).
		- Biological benefit: Prevents one hyper-cellular surgical sample from dominating downstream PCA/UMAP clustering and CellChat interaction probabilities.
3. Cohort Merging & Ribosomal QC (`percent.ribo`) (Steps 4 & 5):
	- Merges all patient Seurat objects into a single unified BCBM object (`bcbm_seurat`), prefixing cell barcodes with sample_id (`add.cell.ids = names(seurat_list)`) to guarantee zero barcode collisions across patients.
	- Computes ribosomal RNA content (percent.ribo via ^RP[SL]) as an additional metabolic/translation QC covariate and saves the uncompressed object (01_filtered_merged_seurat.rds) for fast disk I/O in Phase 2.
