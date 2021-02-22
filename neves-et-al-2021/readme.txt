This folder is self-contained. The R scripts here are enough to reproduce the figures in the Neves et al., 2021 article, based on the model in this repo.

To generate the figures, run the script plot_figures.R, indicating the location of the data in the beginning of the script.

The data is available at: https://osf.io/z2hn9/files/, under the path /Article - Neves et al. 2021/Data

Here's the output of R sessionInfo() when the figures were generated:

R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=pt_BR.UTF-8      
 [2] LC_NUMERIC=C              
 [3] LC_TIME=pt_BR.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=pt_BR.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=pt_BR.UTF-8      
 [8] LC_NAME=C                 
 [9] LC_ADDRESS=C              
[10] LC_TELEPHONE=C            
[11] LC_MEASUREMENT=pt_BR.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils    
[5] datasets  methods   base     

other attached packages:
[1] openxlsx_4.1.5 glue_1.4.1     purrr_0.3.4   
[4] cowplot_1.0.0 

loaded via a namespace (and not attached):
 [1] zip_2.0.4         Rcpp_1.0.4.6     
 [3] pillar_1.4.4      compiler_3.6.3   
 [5] cellranger_1.1.0  forcats_0.5.0    
 [7] tools_3.6.3       lifecycle_0.2.0  
 [9] tibble_3.0.1      gtable_0.3.0     
[11] pkgconfig_2.0.3   rlang_0.4.6      
[13] rstudioapi_0.11   curl_4.3         
[15] yaml_2.2.1        haven_2.3.1      
[17] rio_0.5.16        dplyr_1.0.0      
[19] generics_0.0.2    vctrs_0.3.1      
[21] hms_0.5.3         tidyselect_1.1.0 
[23] grid_3.6.3        data.table_1.13.0
[25] R6_2.4.1          readxl_1.3.1     
[27] foreign_0.8-76    carData_3.0-4    
[29] ggplot2_3.2.1     car_3.0-8        
[31] magrittr_1.5      scales_1.1.1     
[33] ellipsis_0.3.1    abind_1.4-5      
[35] colorspace_1.4-1  stringi_1.4.6    
[37] lazyeval_0.2.2    munsell_0.5.0    
[39] crayon_1.3.4  