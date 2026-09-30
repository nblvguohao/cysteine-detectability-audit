# External inputs fetched for the 2026-09-30 revision analyses

All files below were downloaded on 2026-09-30 into `revision_2026-09-30/external/`; they are not redistributed here and can be re-fetched from the URLs.
Downloads were approved by the authors before they were made.

| file | source URL | bytes | sha256 | note |
|---|---|---|---|---|
| UP000000589_10090.fasta.gz | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz | 8053531 | c24568f70ce9ffd7693a1ecd31c23eaa97a1f43d9305e2ff5401bdf30e2a0a93 | UniProt reference proteomes release 2026_03 (mouse). The original analyses used a 2025 download of the same proteome (sha256 033ac9c4...); identified peptides are checked against the new sequences at their stated positions. |
| UP000059680_39947.fasta.gz | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000059680/UP000059680_39947.fasta.gz | 9338783 | 0362824104575e680053e3168a85f72d0cd5e1f1c96100753510a06ba5f2866b | UniProt 2026_03, rice (Oryza sativa subsp. japonica); covers ~45% of the deposit's leading razor accessions, the rest are fetched per accession. |
| 20240326_061934_xyj_proteome_Report.tsv | ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD072035/20240326_061934_xyj_proteome_Report.tsv | 7765132 | 36ce465600c72a96763aa7522fd7d3d019fed0c16e0de70f396437e508bf955a | PRIDE PXD072035 (rice total proteome, Spectronaut protein-group report) |
| SS-all-peptides.tsv | ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD072089/SS-all-peptides.tsv | 2625176 | 57c4ce64b5747eb591e879565916c04e0901ea228c20e511b2e1804fb384890c | PRIDE PXD072089 (rice persulfidome, MaxQuant peptides table) |
| pLMSNOSite/ (git clone) | https://github.com/KCLabMTU/pLMSNOSite at commit e9158af | ~252 MB working tree | blob sha256: data/train/sequence_train.csv a3b4a1fe...bd0f, data/test/sequence_test.csv d4498e6f...5551 | Blob hashes match Supplemental Note 5. The Windows checkout converts line endings, so hash the git blob (`git show HEAD:path`), not the working-tree file. |
| .venv_tf/ | PyPI | ~0.5 GB | n/a | tensorflow==2.15.1, numpy<2, pandas, scikit-learn, biopython, h5py; used only to run the released pLMSNOSite models |
