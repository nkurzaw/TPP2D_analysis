Analysis of GTP dataset
================
25 June, 2020

# Step-by-step walk through the anlysis

``` r
# This script uses the development version of TPP2D
if(require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("nkurzaw/TPP2D")
```

Load required libraries

``` r
library(TPP2D)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(UpSetR)
```

Define plot style

``` r
theme_paper <- theme_bw(base_size = 6) +
  theme(legend.background = element_blank(), 
        legend.key = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "grey92", size = 0.25),
        panel.grid.minor = element_line(colour = "grey92", size = 0.15),
        panel.border = element_blank(), 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        complete = TRUE,
        axis.line = element_line(color = "black", size = 0.25),
        text = element_text(size = 7),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
```

Annotate GTP and ATP binding proteins

``` r
all_atp_binder <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, 
  keys = "GO:0005524", 
  columns = c("SYMBOL", "IPI"), 
  keytype = "GOALL")
```

    ## 

    ## 'select()' returned 1:many mapping between keys and columns

``` r
all_gtp_binder <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, 
  keys = "GO:0005525", 
  columns = c("SYMBOL", "IPI"), 
  keytype = "GOALL")
```

    ## 'select()' returned 1:many mapping between keys and columns

Download the supplementary table from the journal’s website

``` r
if(!file.exists("Supplementary_Data_3.xlsx")){
  download.file(
    url = "https://www.biorxiv.org/content/biorxiv/early/2020/05/09/2020.05.08.083709/DC2/embed/media-3.xlsx?download=true",
    destfile = "Supplementary_Data_3.xlsx",
    mode = "wb")
}
```

Read in the data and reformat to a data frame as would be obtained after
import of the raw data:

``` r
gtp_raw <- read_xlsx("Supplementary_Data_3.xlsx", sheet = "GTP") %>% 
  dplyr::select(representative,
                clustername,
                qupm,
                qusm,
                experiment,
                temperature,
                matches("sumionarea"),
                matches("rel_fc_protein"))  %>% 
  gather(key, value, matches("sumionarea"), matches("rel_fc_protein")) %>% 
  mutate(conc = as.numeric(gsub("uM", "", gsub(".+_protein_[0-9,H,L]+_[0-9,H,L]+_", "", key))),
         temperature = as.numeric(gsub("C", "", temperature)),
         key = case_when(grepl("sumionarea", key) ~ "raw_value",
                         grepl("rel_fc", key) ~ "rel_value")) %>% 
  spread(key, value) %>% 
  arrange(representative, temperature, conc) %>% 
  group_by(clustername, temperature, conc) %>% 
  filter(qupm == max(qupm), 
         qusm == max(qusm), 
         raw_value == max(raw_value)) %>% 
  filter(!duplicated(clustername)) %>% 
  ungroup %>% 
  mutate(log2_value = log2(raw_value),
         log_conc = log10(conc/1e6)) %>% 
  filter(qupm > 1)

# resolve ambiguous protein names
gtp_fil <- resolveAmbiguousProteinNames(gtp_raw)
  
# recompute reporter ion signal from robust Isobarquant fold changes
gtp_df <- recomputeSignalFromRatios(gtp_fil)
```

Compute null and alternative model fits and extract parameters

``` r
gtp_params_df <- getModelParamsDf(gtp_df, maxit = 500)
saveRDS(gtp_params_df, file = "../pre_run_data/gtp_params_df.rds")
```

Compute *F* statistics

``` r
gtp_fstat_df <- computeFStatFromParams(gtp_params_df)
```

Get \(B\) datasets expected under the null model and perform model
fitting and compute F statistics to obtain a null distribution for FDR
calibration:

``` r
set.seed(12, kind = "L'Ecuyer-CMRG")
gtp_null_df <- bootstrapNullAlternativeModel(
  df = gtp_df, params_df = gtp_params_df, 
  maxit = 500, B = 100,
  BPPARAM = BiocParallel::MulticoreParam(workers = 20, progressbar = TRUE),
  verbose = FALSE)
saveRDS(jq1_null_df, file = "../pre_run_data/jq1_null_df.rds")
```

Compute FDR and find hits:

``` r
gtp_fdr_df <- getFDR(df_out = gtp_fstat_df,
                     df_null = gtp_null_df,
                     squeezeDenominator = TRUE)
  
gtp_hits_df <- findHits(gtp_fdr_df, alpha = 0.1)
```

``` r
ggplot(gtp_fdr_df %>% 
         filter(dataset == "true"), 
       aes(log2(rssH0 - rssH1), asinh(F_statistic))) +
  geom_point(color = "gray", alpha = 0.5, size = 0.5) + 
  geom_point(color = "black", alpha = 0.5, 
             size = 0.5,
             data = gtp_hits_df %>% 
              filter(!clustername %in% 
                       all_atp_binder$SYMBOL,
                     !clustername %in% 
                       all_gtp_binder$SYMBOL)) + 
  geom_point(color = "darkgreen", alpha = 0.5, 
             size = 0.5,
             data = filter(gtp_hits_df, clustername %in% 
                             all_atp_binder$SYMBOL)) + 
  geom_point(color = "violet", alpha = 0.5, 
             size = 0.5,
             data = filter(gtp_hits_df, clustername %in% 
                             all_gtp_binder$SYMBOL)) +
  facet_wrap(~nObsRound, scales = "free", ncol = 5) +
   labs(x = expression('log'[2]~'(RSS'^0~' - RSS'^1~')'),
       y = expression('asinh('*italic(F)*' statistic)')) +
  coord_cartesian(xlim = c(-12.5, 7.5), ylim = c(0, 6)) +
  ggtitle("GTP gelfil. lysate experiment") +
  theme_paper
```

<img src="md_files/gtp-unnamed-chunk-13-1.png" width="100%" />

# Compare sets

``` r
gtp_threshold_df <- read_xlsx("Supplementary_Data_3.xlsx", sheet = "GTP") %>% 
  dplyr::select(clustername, stabilized_hits_found_by_threshold) %>% 
  filter(stabilized_hits_found_by_threshold)

gtp_hits_10per <- findHits(gtp_fdr_df, 0.10) %>% 
  filter(detected_effectH1 == "stability", slopeH1 > 0)
gtp_hits_1per <- findHits(gtp_fdr_df, 0.01) %>% 
  filter(detected_effectH1 == "stability", slopeH1 > 0)


set_intersect_df <- data.frame(
  gene_name = filter(gtp_fdr_df, !duplicated(clustername))$clustername) %>% 
  mutate(`ATP binders` = as.numeric(gene_name %in% all_atp_binder$SYMBOL),
         `GTP binders` = as.numeric(gene_name %in% all_gtp_binder$SYMBOL),
         `DLPTP hits at 10% FDR` = 
           as.numeric(gene_name %in% gtp_hits_10per$clustername),
         `DLPTP hits at 1% FDR` = 
           as.numeric(gene_name %in% gtp_hits_1per$clustername),
         `Hits found by threshold-based approach` = 
           as.numeric(gene_name %in% gtp_threshold_df$clustername))

upset(set_intersect_df, nsets = 5)
```

<img src="md_files/gtp-unnamed-chunk-14-1.png" width="100%" />

``` r
sessionInfo()
```

    ## R version 4.0.0 Patched (2020-05-04 r78358)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] UpSetR_1.4.0  readxl_1.3.1  ggplot2_3.3.2 tidyr_1.1.0   TPP2D_1.5.5  
    ## [6] dplyr_1.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0     xfun_0.14            purrr_0.3.4         
    ##  [4] colorspace_1.4-1     vctrs_0.3.0          generics_0.0.2      
    ##  [7] htmltools_0.4.0      stats4_4.0.0         yaml_2.2.1          
    ## [10] blob_1.2.1           rlang_0.4.6          pillar_1.4.4        
    ## [13] glue_1.4.1           withr_2.2.0          DBI_1.1.0           
    ## [16] BiocParallel_1.22.0  BiocGenerics_0.34.0  bit64_0.9-7         
    ## [19] foreach_1.5.0        lifecycle_0.2.0      plyr_1.8.6          
    ## [22] stringr_1.4.0        munsell_0.5.0        gtable_0.3.0        
    ## [25] cellranger_1.1.0     zip_2.0.4            codetools_0.2-16    
    ## [28] evaluate_0.14        memoise_1.1.0        labeling_0.3        
    ## [31] Biobase_2.48.0       knitr_1.28           IRanges_2.22.2      
    ## [34] doParallel_1.0.15    parallel_4.0.0       AnnotationDbi_1.50.0
    ## [37] Rcpp_1.0.4.6         scales_1.1.1         limma_3.44.1        
    ## [40] org.Hs.eg.db_3.11.4  S4Vectors_0.26.1     farver_2.0.3        
    ## [43] bit_1.1-15.2         gridExtra_2.3        digest_0.6.25       
    ## [46] stringi_1.4.6        openxlsx_4.1.5       grid_4.0.0          
    ## [49] tools_4.0.0          bitops_1.0-6         magrittr_1.5        
    ## [52] RCurl_1.98-1.2       tibble_3.0.1         RSQLite_2.2.0       
    ## [55] crayon_1.3.4         pkgconfig_2.0.3      MASS_7.3-51.6       
    ## [58] ellipsis_0.3.1       rmarkdown_2.2        iterators_1.0.12    
    ## [61] R6_2.4.1             compiler_4.0.0
