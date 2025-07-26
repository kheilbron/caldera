# CALDERA (Calling disease-related genes)

## General overview
CALDERA is a tool for identifying causal genes in GWAS loci. Many alternative tools (*e.g.*, L2G, Ei, FLAMES) use complex XGBoost models with >45 predictive features. CALDERA achieves similar or better performance despite using a simple logistic regression model with only 4 features and 3 covariates. This makes it much easier to understand why CALDERA chooses the genes that it does. Unlike other methods, CALDERA attempts for potential biases in its ground truth training data. Full details can be found in the [preprint](https://www.medrxiv.org/content/10.1101/2024.07.26.24311057v2).

This repo contains everything needed to run CALDERA given: 1) a PoPS output file, and 2) a file containing credible set information.

## How to run CALDERA
First, clone this repo
```
git clone https://github.com/kheilbron/caldera.git
```

Second, open an R session and source the z_caldera.R script
```
caldera_path <- "your_local_path_to_the_caldera_repo"
source( file.path( caldera_path, "z_caldera.R" ) )
```

Finally, run CALDERA by providing the two required input files (described in more detail in the next section)
```
pops_file  <- "path_to_your_PoPS_output_file"
cs_file    <- "path_to_your_credible_set_file"
caldera( pops_file    = pops_file,
         cs_file      = cs_file,
         caldera_path = caldera_path )
```

## CALDERA input files
| File       | Description |
| ---------- | ---------- |
| pops_file  | PoPS output file (pops.preds). Must contain columns: ENSGID and PoPS_Score. |
| cs_file    | A (comma/tab/etc. delimited) file containing all credible set (CS) variants across all loci. Used to define loci, the genes therein, and coding PIPs for these genes. One CS variant per row. Must contain columns: "chr" (chromosome), "bp" (GRCh37 position), and "pip" (posterior inclusion probability). Must also contain a "locus" column identifying which credible set each variant belongs to (can be an integer or a string). Optionally, you can include a "c_gene" column containing the ENSGID of the gene whose coding sequence is affected by the given variant (and NA if the variant does not affect any coding sequences). Alternatively, an "rsid" column containing rsIDs can be provided and variants will automatically be annotated using a list of ~76k UK Biobank coding variants. If "rsid" is not provided, variants will be automatically annotated based on "chr" and (GRCh37) "bp". |

## CALDERA output
The caldera() R function returns a data.table containing the following columns:

| Column        | Description |
| ----------    | ---------- |
| locus         | The locus identifier, taken from the input credible set file |
| locus_pos     | A string containing the chromosome, left boundary, and right boundary of the locus (each separated by underscores) |
| gene          | Gene name |
| caldera       | CALDERA-predicted probability that this gene is causal for the GWAS trait |
| multi         | As above, but without normalizing probabilities to sum to 100% in each locus |
| n_genes       | Number of protein-coding genes in the locus |
| dist          | Distance between the lead GWAS variant and any part of the gene (bp) |
| pops          | PoPS score for the gene |
| coding        | Sum of posterior inclusion probabilities of coding variants in this credible set affecting this gene |
| ensgid        | ENSEMBL gene ID |


