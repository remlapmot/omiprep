# Export Data

`Metaboprep` can export data to various formats.

## Setup

### load the metaboprep library

``` r
library(metaboprep)
```

### Read in the data and make a Metaboprep object

Create a `Metaboprep` object as described in the [Getting
Started](https://mrcieu.github.io/metaboprep/articles/index.md)
vignette.

``` r
# read in the metabolon data as a list object
datain     <- read_metabolon(system.file("extdata", "metabolon_v1.1_example.xlsx", package = "metaboprep"), 
                            sheet="OrigScale", 
                            return_Metaboprep = FALSE)

# build the Metaboprep class object
mydata      <- Metaboprep(data = datain$data, samples = datain$samples, features = datain$features)
```

## Run the quality control

``` r
## Adding suppressWarnings() to avoid deparse() error when rendering vignette with S7 method warnings
mydata         <- suppressWarnings( quality_control(mydata, cores = 1) )
#> 
#> ── Starting Metabolite QC Process ──────────────────────────────────────────────
#> ℹ Validating input parameters
#> ✔ Validating input parameters [8ms]
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> AF =  2
#> ✔ Sample & Feature Summary Statistics for raw data [579ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [35ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [21ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [17m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 2 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> 
#> ℹ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> ✔ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [26ms]
#> 
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> AF =  2
#> ! The stated max PCs [max_num_pcs=10] to use in PCA outlier assessment is greater than the number of available informative PCs [2]
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…✔ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> 
#> ℹ Creating final QC dataset...
#> AF =  2
#> ✔ Creating final QC dataset... [634ms]
#> 
#> ℹ Metabolite QC Process Completed
#> 
#> ℹ Metabolite QC Process Completed── Step timings ──
#> ℹ Metabolite QC Process Completed
#> ℹ Metabolite QC Process Completed
#>                         step seconds  pct
#>                   validation    0.00  0.0
#>                summarise_raw    0.57 27.8
#>                   copy_layer    0.00  0.0
#>   extreme_sample_missingness    0.00  0.0
#>  extreme_feature_missingness    0.00  0.0
#>           sample_missingness    0.00  0.0
#>          feature_missingness    0.00  0.0
#>              total_peak_area    0.00  0.0
#>                summarise_pca    0.65 31.7
#>              summarise_final    0.60 29.2
#>                        total    2.05 99.9
#> ✔ Metabolite QC Process Completed [40ms]
```

## Export Metaboprep

``` r
# where to put the files
output_dir <- file.path(getwd(), "output")

# run export
export(mydata, directory = output_dir, format = "metaboprep")
#> Exporting in metaboprep format to: 
#>      /home/runner/work/metaboprep/metaboprep/vignettes/output

# view output directory files
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
unname(sapply(files, function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  paste(tail(parts, 4), collapse = .Platform$file.sep)
}))
#>  [1] "output/metaboprep_export_2026_04_17/input/config.yml"         
#>  [2] "output/metaboprep_export_2026_04_17/input/data.tsv"           
#>  [3] "output/metaboprep_export_2026_04_17/input/feature_summary.tsv"
#>  [4] "output/metaboprep_export_2026_04_17/input/features.tsv"       
#>  [5] "output/metaboprep_export_2026_04_17/input/sample_summary.tsv" 
#>  [6] "output/metaboprep_export_2026_04_17/input/samples.tsv"        
#>  [7] "output/metaboprep_export_2026_04_17/qc/config.yml"            
#>  [8] "output/metaboprep_export_2026_04_17/qc/data.tsv"              
#>  [9] "output/metaboprep_export_2026_04_17/qc/feature_summary.tsv"   
#> [10] "output/metaboprep_export_2026_04_17/qc/feature_tree.RDS"      
#> [11] "output/metaboprep_export_2026_04_17/qc/features.tsv"          
#> [12] "output/metaboprep_export_2026_04_17/qc/sample_summary.tsv"    
#> [13] "output/metaboprep_export_2026_04_17/qc/samples.tsv"           
#> [14] "output/metaboprep_export_2026_04_17/qc/var_exp.tsv"
```
