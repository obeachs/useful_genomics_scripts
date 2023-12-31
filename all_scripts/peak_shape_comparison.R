library(rtracklayer)
library(parallel)
library(BH)
library(Rcpp)
library(FunChIP)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(zoo)
data(GR100)
library(ggpmisc)
library(DiffBind)
library(pracma)
library(rstatix)
library(ggpubr)
source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')


load("~/SIC-ChIP/SIC-ChIP_functions.RData") 

clf28_bam <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/CHKPEI85219060189-6_190716_X602_FCH2TJYCCX2_L4_CHKPEI85219060189-6_phix_removed.bam'
clf28_alp1_bam <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/CHKPEI85219060188-11_190716_X602_FCH2TJYCCX2_L3_CHKPEI85219060188-11_phix_removed.bam'
clf28_alp2_bam <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/CHKPEI85219060188-6_190716_X602_FCH2TJYCCX2_L3_CHKPEI85219060188-6_phix_removed.bam'
clf_bw <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw')
double_peaks <- list('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/all_rescues/clf_alp_rescue_peaks/col-0_peaks_rescue_genes.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/all_rescues/clf_alp_rescue_peaks/clf28_peaks_rescue_genes.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/all_rescues/clf_alp_rescue_peaks/clf28_alp1_peaks_rescue_genes.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/all_rescues/clf_alp_rescue_peaks/clf28_alp2_peaks_rescue_genes.csv')
double_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_H3K27me3_RPGC.bw')
true_peaks <- list('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/col-0_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp1_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp2_sig_peaks_rescues.csv')
all_peaks <- list('/Volumes/sesame/ALP_Omics/ChIP/validations/col_0_chr5_H3K27me3_sig_peaks.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf28_chr5_H3K27me3_sig_peaks.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf28_alp1_chr5_H3K27me3_sig_peaks.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf28_alp2_chr5_H3K27me3_sig_peaks.csv')


count_local_peaks <- function(local_maxima){
            #Order the dataframe by the start positions
            local_maxima <- local_maxima[order(local_maxima$start), ]
            # Initialize counters
            count <- 0
            current_max_end <- -Inf

            # Loop through the dataframe and compare the start with 
            # the current 'max_end' variable. Once the start is greater than
            # the next end (as soon as it starts to curve downwards) we take
            # that as a local maximum. The 500 is to make sure that there are 
            # no two local maxima called within 500bp of each other.
            for (i in 1:nrow(local_maxima)) {
            if (local_maxima$start[i] - current_max_end > 500) {
                # This is a new local maximum
                count <- count + 1
            }
            # Update the current maximum's end position
            current_max_end <- max(current_max_end, local_maxima$end[i])
            }
            count
    }


count_local_peaks_2 <- function(local_maxima){}

extract_peak_region <- function(bigwig_df, chrom, start, end){
    bigwig_df$seqnames <- gsub('Chr','',bigwig_df$seqnames)
    chrom_name <- chrom
    subset_df <- bigwig_df[bigwig_df$seqnames == chrom_name,]
    #start <- bigwig_df$start[1]
    #end <- bigwig_df$end[1]
    #This next two lines are finding the start and end 
    # coordinates that are as close as possible to the 
    # start and end denoted by the peak region - they often 
    # differ by a little bit because of the different window
    # sizes used for each
    seq_start <- which(abs(subset_df$start - start) == min(abs(subset_df$start - start)))
    seq_end <- which(abs(subset_df$start - end) == min(abs(subset_df$start - end)))
    subset_df <- subset_df[seq_start:seq_end,]
    subset_df
}


smooth_peak_region <- function(peak_bw_df){
    #Not sure if this is necesary but convering it to granges object
    # before the transformation
    subsetted_gr <- with(peak_bw_df, GRanges(seqnames, IRanges(start, end), score = score))
    #Setting the rolling window for thing to be 200, seems to be in line
    # with other peak calling methods
    thed_gr <- as.data.frame(transform(subsetted_gr, score_thed = rollmean(score, k = 100, align = "center", fill = NA)))
    thed_gr[is.na(thed_gr)] <- 0
    thed_gr
}

extract_local_maxima <- function(raw_peaks,th_peaks){
    diff1 <- c(NA, diff(th_peaks$score_thed))
    diff2 <- c(diff(th_peaks$score_thed), NA)
    is_local_max <- diff1 > 0 & diff2 < 0
    is_local_max[is.na(is_local_max)] <- FALSE
    local_maxima <- raw_peaks[is_local_max,]
    local_maxima
}

extract_local_maxima_raw_only <- function(raw_peaks){
    diff1 <- c(NA, diff(raw_peaks$score))
    diff2 <- c(diff(raw_peaks$score), NA)
    is_local_max <- diff1 > 0 & diff2 < 0
    is_local_max[is.na(is_local_max)] <- FALSE
    local_maxima <- raw_peaks[is_local_max,]
    local_maxima
}

get_peak_info <- function(peaks, bigwig, outname){
    results <- data.frame(
    peak_id = integer(0),
    max_height = numeric(0),
    num_peaks = integer(0),
    area = numeric(0),
    width = integer(0)
    )


    peakfile <- unique(read.csv(peaks))
    print(head(peakfile))
    bw <- plotgardener::readBigwig(bigwig)

    for(i in 1:length(peakfile$chr)){
        # Code that might generate an error
        print(peakfile$feature[i])
        print(i)
        p_start <- peakfile$start[i]
        p_end <- peakfile$end[i]
        p_chr <- peakfile$chr[i]
        p_chr <- gsub('Chr','', p_chr)
        print(p_start)
        print(p_end)
        print(p_chr)
        print(peakfile$width[i])
        if (as.numeric(peakfile$width[i]) > 1000){
            print('extract_peak_region')
            peak_bw <- extract_peak_region(bigwig_df=bw, start = p_start, end=p_end, chrom = p_chr)
            print('thing peaks')
            smooth_peak_bw <- smooth_peak_region(peak_bw)
            print(min(smooth_peak_bw$start))
            print(max(smooth_peak_bw$end))
            print('getting local maxima')
            local_peaks <- extract_local_maxima(peak_bw, smooth_peak_bw)
            print(local_peaks)
            
            peak_name <- peakfile$peak[i]
            #local_count <- count_local_peaks(local_peaks)
            local_count <- tryCatch(
            {
            count_local_peaks(local_peaks)
            },
            error = function(e) {
            # Handle the error by returning 1
            return(1)
            }
        )
            area <- sum(peak_bw$score)
            max_h <- max(smooth_peak_bw$score_thed)
            width <- peakfile$width[i]
            w_h_ratio <- (width/max_h)
            output_list <- c(peak_name, max_h, local_count, area, width, w_h_ratio)
            results <- rbind(results, output_list)
        }
        else{
            peak_name <- peakfile$peak[i]
            local_count <- 1
            smooth_peak_bw <- smooth_peak_region(peak_bw)
            area <- sum(peak_bw$score)
            max_h <- max(smooth_peak_bw$score_thed)
            width <- peakfile$width[i]
            w_h_ratio <- (width/max_h)
            output_list <- c(peak_name, max_h, local_count, area, width,w_h_ratio)

            results <- rbind(results, output_list)

        }
    }
    names(results) <- c('peak_id','max_height','num_peaks','area','width', 'w_h_ratio')
    write.csv(results, outname)
}



get_peak_info_raw_peaks <- function(peaks, bigwig, outname){

    results <- data.frame(
    peak_id = integer(0),
    max_height = numeric(0),
    num_peaks = integer(0),
    area = numeric(0),
    width = integer(0)
    )
    peakfile <- read.csv(peaks)
    peakfile <- mutate(peakfile,peak=rownames(peakfile))
    print(head(peakfile))
    bw <- plotgardener::readBigwig(bigwig)
    for(i in 1:length(peakfile$chr)){
        # Code that might generate an error
        p_start <- peakfile$start[i]
        p_end <- peakfile$end[i]
        p_chr <- peakfile$chr[i]
        p_chr <- gsub('Chr','', p_chr)
        print(peakfile$width[i])
        print(p_start)
        print(p_end)
        if (as.numeric(peakfile$width[i]) > 1000){
            print('extract_peak_region')
             peak_bw <- extract_peak_region(bigwig_df = bw, start = p_start, end = p_end, chrom = p_chr)
            print('thing peaks')
            smooth_peak_bw <- smooth_peak_region(peak_bw)
            print(min(smooth_peak_bw$start))
            print(max(smooth_peak_bw$end))
            print('getting local maxima')
            local_peaks <- extract_local_maxima(peak_bw, smooth_peak_bw)
            print(local_peaks)
            
            peak_name <- peakfile$peak[i]
            #local_count <- count_local_peaks(local_peaks)
            local_count <- tryCatch(
            {
            count_local_peaks(local_peaks)
            },
            error = function(e) {
            # Handle the error by returning 1
            return(1)
            }
        )
            area <- sum(peak_bw$score)
            max_h <- max(smooth_peak_bw$score_thed)
            width <- peakfile$width[i]
            ratio <- width/max_h
            output_list <- c(peak_name, max_h, local_count, area, width,ratio)
            results <- rbind(results, output_list)
        }
        else{
            peak_bw <- extract_peak_region(bigwig_df=bw, start = p_start, end=p_end, chrom = p_chr)
            smooth_peak_bw <- smooth_peak_region(peak_bw)
            peak_name <- peakfile$peak[i]
            local_count <- 1
            area <- sum(peak_bw$score)
            max_h <- max(smooth_peak_bw$score_thed)
            width <- peakfile$width[i]
            ratio <- width/max_h
            output_list <- c(peak_name, max_h, local_count, area, width, ratio)
            results <- rbind(results, output_list)

        }
    }
    names(results) <- c('peak_id','max_height','num_peaks','area','width', 'w_h_ratio')
    write.csv(results, outname)
}


for(i in length(double_bigwigs)){
    print(all_peaks[[i]])
    out <- paste0(tools::file_path_sans_ext(all_peaks[[4]]), '_stats')
    get_peak_info_raw_peaks(peaks=all_peaks[[4]],bigwig=double_bigwigs[[4]], outname = out)
}


out <- paste0(tools::file_path_sans_ext(all_peaks[[4]]), '_stats')
get_peak_info_raw_peaks(peaks=all_peaks[[4]],bigwig=double_bigwigs[[4]], outname = out)



col <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/', pattern='col_0', full.names = T)
col_l <- lapply(col, read.csv)
col <- dplyr::bind_rows(col_l) %>%mutate(exp='col-0')
clf <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/', pattern='clf28_', full.names = T)
clf_l <- lapply(clf, read.csv)
clf <- dplyr::bind_rows(clf_l) %>%mutate(exp='clf28')
alp1 <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/', pattern='alp1', full.names = T)
alp1_l <- lapply(alp1, read.csv)
alp1 <- dplyr::bind_rows(alp1_l) %>%mutate(exp='clf28_alp1')
alp2 <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/', pattern='alp2', full.names = T)
alp2_l <- lapply(alp2,read.csv)
alp2 <- dplyr::bind_rows(alp2_l) %>%mutate(exp='clf28_alp2')


col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/stats/col-0_sig_peaks_rescues_stats_rpgc') %>% mutate(exp='col-0')
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/stats/clf28_sig_peaks_rescues_stats_rpgc')%>% mutate(exp='clf28')
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/stats/clf28_alp1_sig_peaks_rescues_stats_rpgc')%>% mutate(exp='clf28_alp1')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/stats/clf28_alp2_sig_peaks_rescues_stats_rpgc')%>% mutate(exp='clf28_alp2')









df <- rbind(col,clf,alp1,alp2)
kmeans(as.matrix())
pwc <- df %>%
  pairwise_t_test(area ~ exp, p.adjust.method = "bonferroni") %>% rename(area=.y.)

p <- ggplot(df, aes(x=exp, y=width, fill=exp)) + 
  geom_violin() + scale_fill_brewer(palette="Paired") + theme + theme(legend.text=element_text(size=20))+
    stat_summary(fun.y=median, geom="point", size=4, color="red")+
    stat_summary(fun.y=mean, geom="point", size=4, color='white')+ theme(axis.text=element_text(size=15),
        axis.title=element_text(size=14,face="bold"))
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/stats/chr2_plots/width.png',p)



my_comparisons <- list( c("clf28", "col-0"), c("clf28", "clf28_alp1"), c("clf28", "clf28_alp2"),
c("clf28_alp1", "clf28_alp2"),c("col-0", "clf28_alp2"),c("col-0", "clf28_alp1"))

d <- ggboxplot(df, x = "exp", y = "w_h_ratio",color = "exp", palette = "Paired")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1)+ theme(axis.text=element_text(size=15),
        axis.title=element_text(size=14,face="bold"))


d <- ggviolin(
  df, x = "exp", y = "width", color = "exp", palette = "Paired",
  add = "jitter", # Add jittered points for better visualization
  fill = "exp"    # Fill the violins with colors based on the "exp" variable
    ) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_summary(
    fun = "mean",   # Calculate mean for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of mean points
    color = "black" # Set the color of mean points
  ) +
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/stats/cpm/max_height.png',d)



source('/Volumes/sesame/ALP_Omics/ChIP/scripts/annotate_peaks_transposable_elements.R')
annotate_peaks_with_transposons('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/col-0_peaks.broadPeak',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/col-0_peaks_annotated','Chromosome','Start','End')
###Looking at the ref6-c alp2 double mutants (relative to ref6)
peaks <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/peaks', pattern='control_peaks.csv', full.names=T)
bigiwigs <-c('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/col_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_rpkm.bw')


macs_col <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/col-0_peaks_annotated.csv')
macs_clf <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/clf28_peaks_annotated.csv')
macs_clf_ <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/clf28_alp2_peaks_annotated.csv')

refpeaks <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% filter(FC_KO > 1.5 | FC_KO < 0.75)
write.csv(refpeaks,'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_alp2_ref6_control_sig_peaks', quote=F, row.names=F)

for(i in 1:length(peaks)){
    print(peaks[[i]])
    out <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control//true_peaks/stats/',basename(peaks[[i]]), '_stats')
    print(out)
    get_peak_info_raw_peaks(peaks='/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/true_peaks.csv',bigwig=bigiwigs[[i]], outname = out)
}
get_peak_info_raw_peaks(peaks='/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_sig_peaks.csv',
bigwig=bigiwigs[[4]], outname = '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_sig_peaks_stats')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_annotated_midpeak_transposons_fixed.csv') %>%
filter(feature %in% refpeaks$feature)
print(length(alp2$feature))

raw_peaklist <- list('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/col-0_me3_annotated_midpeak_transposons_fixed.csv',
'/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_annotated_midpeak_transposons_fixed.csv',
'/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_annotated_midpeak_transposons_fixed.csv',
'/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_annotated_midpeak_transposons_fixed.csv')
for(i in raw_peaklist){
    print(tools::file_path_sans_ext(i))
    peak <- read.csv(i) %>% filter(feature %in% refpeaks$feature)
    write.csv(peak, paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/',basename(tools::file_path_sans_ext(i)),'_ref6_control_peaks.csv'), quote=F, row.names=F)
}

col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/col-0_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='Col-0')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/alp2_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='alp2-1')
ref6alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/ref6_alp2_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='ref6_alp2-1')
ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/ref6_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='ref6')

df <- rbind(col,alp2,ref6alp2,ref6) %>% mutate(width=ifelse(max_height==0,0,width))
my_comparisons <- list( c("Col-0", "alp2-1"), c("Col-0", "ref6_alp2-1"), c("Col-0", "ref6"),
c("ref6", "alp2-1"),c("ref6_alp2-1", "alp2-1"),c("ref6_alp2-1", "ref6"))

my_comparisons_no_col <- list(c("ref6", "alp2-1"),c("ref6_alp2-1", "alp2-1"),c("ref6_alp2-1", "ref6"))

df

d <- ggviolin(
  df, x = "exp", y = "width", color = "exp", palette = "Paired",
  add = "jitter", # Add jittered points for better visualization
  fill = "exp"    # Fill the violins with colors based on the "exp" variable
    ) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_summary(
    fun = "mean",   # Calculate mean for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of mean points
    color = "black" # Set the color of mean points
  ) +
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/plots/width.png',d)







###alp2_single mutant checks
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/col-0_me3_annotated_midpeak.csv') %>% filter(log2FoldChange >2)%>% filter(FDR < 0.05)
write.csv(col, '/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/col_raw_sig_peaks.csv', quote=F, row.names=F)
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_annotated_midpeak.csv')%>% filter(log2FoldChange > 2)%>% filter(FDR < 0.05)
write.csv(col, '/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp1_raw_sig_peaks.csv', quote=F, row.names=F)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_annotated_midpeak.csv')%>% filter(log2FoldChange >2) %>% filter(FDR < 0.05)
write.csv(col, '/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp2_raw_sig_peaks.csv', quote=F, row.names=F)

peaks <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/', pattern='sig_peaks.csv',full.names = T)
bigwigs <- c('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp1_rpkm.bw', '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/col_rpkm.bw')
for(i in 1:length(peaks)){
    print(peaks[[i]])
    out <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/',basename(peaks[[i]]), '_stats')
    print(out)
    get_peak_info_raw_peaks(peaks=peaks[[i]],bigwig=bigwigs[[i]], outname = out)
}


##Getting the features where alp1 and alp2 are significantly different to col
alp1_dif <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO > 1.5)
alp2_dif <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1.5)
genes <- c(alp1_dif$feature, alp2_dif$feature)
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/col-0_me3_annotated_midpeak.csv') %>% filter(feature %in% genes)
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_annotated_midpeak.csv')%>% filter(log2FoldChange > 2)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_annotated_midpeak.csv')%>% filter(feature %in% genes)
diffpeaks <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks', full.names = T)
for(i in 1:length(diffpeaks)){
    print(diffpeaks[[i]])
    out <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/stats',basename(peaks[[i]]), '_stats')
    print(out)
    get_peak_info_raw_peaks(peaks=diffpeaks[[i]],bigwig=bigwigs[[i]], outname = out)
}




col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/stats/statscol_raw_sig_peaks.csv_stats') %>% mutate(exp='Col-0')
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/stats/statsalp1_raw_sig_peaks.csv_stats') %>% mutate(exp='alp1-1')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/stats/statsalp2_raw_sig_peaks.csv_stats') %>% mutate(exp='alp2-1')

df <- rbind(col, alp1, alp2)
my_comparisons <- list( c("Col-0", "alp1-1"), c("Col-0", "alp2-1"), c("alp2-1", "alp1-1"))




d <- ggviolin(
  df, x = "exp", y = "w_h_ratio", color = "exp", palette = "Paired",
  add = "jitter", # Add jittered points for better visualization
  fill = "exp"    # Fill the violins with colors based on the "exp" variable
    ) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_summary(
    fun = "mean",   # Calculate mean for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of mean points
    color = "black" # Set the color of mean points
  ) +
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/stats/plots/ratio.png',d)








source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')
col <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/col_rpkm.bw'
alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2_rpkm.bw'
ref6_alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_alp2_rpkm.bw'
ref6 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_rpkm.bw'

col <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_rpkm.bw'
clf <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw'
clfalp1 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw'
clfalp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_rpkm.bw'
list_2019 <- list(col, clf, clfalp1, clfalp2)
peak_swn_2019 <- list('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/col_peaks_swn_down.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_peaks_swn_down.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_alp1_peaks_swn_down.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_alp2_peaks_swn_down.csv')


col_raw_swn <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/col-0_me3_annotated_midpeak.csv') %>% filter(feature %in% shu_peaks$swn_K27.decrease)
write.csv(col_raw_swn,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/col_peaks_swn_down.csv', quote=F, row.names=F)
clf_raw_swn <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_annotated_midpeak.csv') %>% filter(feature %in% shu_peaks$swn_K27.decrease)
write.csv(clf_raw_swn,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_peaks_swn_down.csv', quote=F, row.names=F)
clfalp1_raw_swn <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_annotated_midpeak.csv') %>% filter(feature %in% shu_peaks$swn_K27.decrease)
write.csv(clfalp1_raw_swn,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_alp1_peaks_swn_down.csv', quote=F, row.names=F)
clfalp2_raw_swn <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_annotated_midpeak.csv') %>% filter(feature %in% shu_peaks$swn_K27.decrease)
write.csv(clfalp2_raw_swn,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/clf_alp2_peaks_swn_down.csv', quote=F, row.names=F)




shu_peaks <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III')
panther_go_maker(clf_peaks_swn_change$feature, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/SWN_decrease_clf_2019_pantherGO')
clf_peaks_swn_change <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_annotated_midpeak.csv') %>% filter(feature %in% shu_peaks$swn_K27.decrease)
names <- read.csv('~/Salba_RNA/genelists/all_hits_names.csv')%>% dplyr::select(tair, gene.name) %>% filter(tair %in% clf_peaks_swn_change$feature) %>%
distinct()
genes_of_interest <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps_with_clf_double/alp2_ref_control_overlap_with_rescues.csv')
#Usage example for plot_peaks
glist <- list(alp2, ref6, ref6_alp2, col)
glist <- list(clf, clfalp1,clfalp2, col)
titles <- list('alp2', 'ref6', 'ref6 alp2','Col-0')
titles <- list('clf28', 'clf28 alp1', 'clf28 alp2','Col-0')



for(i in 1:length(clf_peaks_swn_change$feature)){
  print(i)
  start <- clf_peaks_swn_change$start[i]
  end <- clf_peaks_swn_change$end[i]
  gene <- clf_peaks_swn_change$feature[i]
  print(typeof(start))
  print(typeof(end))
  print(typeof(gene))
  chr <- as.numeric(gsub('Chr','',clf_peaks_swn_change$chr[i]))
  print(typeof(chr))
 plot_peaks(glist,chr,start-2000,end + 2000,gene, titles, paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps_with_clf_double/', gene,'_SWN_decrease.png'))
}

for(i in 1:length(list_2019)){
  print(i)
    out <- paste0(tools::file_path_sans_ext(peak_swn_2019[[i]]), '_stats')
    get_peak_info_raw_peaks(peaks=peak_swn_2019[[i]],bigwig=list_2019[[i]], outname = out)
}

col_swn_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/stats/col_peaks_swn_down_stats') %>% mutate(exp='Col-0')
clf_swn_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/stats/clf_peaks_swn_down_stats') %>% mutate(exp='clf28')
clf28alp1_swn_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/stats/clf_alp1_peaks_swn_down_stats') %>% mutate(exp='clf28 alp1-1')
clf28alp2_swn_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/stats/clf_alp2_peaks_swn_down_stats') %>% mutate(exp='clf28 alp2-1')

df <- rbind(col_swn_stats, clf_swn_stats, clf28alp1_swn_stats,clf28alp2_swn_stats)
my_comparisons <- list( c("Col-0", "clf28 alp1-1"), c("Col-0", "clf28 alp2-1"), c("Col-0", "clf28"), c('clf28 alp1-1','clf28 alp2-1'),c('clf28','clf28 alp1-1'),
c('clf28','clf28 alp2-1'))
sample_sizes <- df %>%
  group_by(exp) %>%
  summarise(sample_size = n())

df <- merge(df, sample_sizes, by = "exp")

d <- ggviolin(
  df, x = "exp", y = "width", color = "exp", palette = "Paired",
  fill = "exp"
    ) +
  stat_summary(
    fun = "mean",   # Calculate mean for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of mean points
    color = "black" # Set the color of mean points
  ) +
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) + theme+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold")
  ) + 
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  )

d_neg <- ggviolin(df, x = "exp", y = "width", palette = "Paired", fill = "exp") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    limits = c(0,max(df$width) + 2000)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    geom = "pointrange",  # Use pointrange for comparison bars
    position = position_dodge(width = 0.75, height = 10000)  # Adjust the width as needed,
  )+
  # Adding sample size label for each 'exp' group
  geom_text(aes(x = exp, y = 0, label = paste("N =", sample_size)), vjust = 2, hjust = 0.5, size = 4)


  ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/stats/width.pdf',
         d_neg, height=10, width=10,limitsize = FALSE)

alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_all_binding_and_down/2020_stats_swn_all/alp1_me3_annotated_midpeak_shared_swn_shared_stats') %>% mutate(exp='alp1')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_all_binding_and_down/2020_stats_swn_all/alp2_me3_annotated_midpeak_shared_swn_shared_stats') %>% mutate(exp='alp2')
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_all_binding_and_down/2020_stats_swn_all/col-0_me3_annotated_midpeak_shared_swn_shared_stats') %>% mutate(exp='Col-0')
my_comparisons <- list(c('alp1','alp2'),c('alp1','Col-0'),c('Col-0','alp2'))
df <- rbind(col, alp1, alp2)
sample_sizes <- df %>%
  group_by(exp) %>%
  summarise(sample_size = n())

df <- merge(df, sample_sizes, by = "exp") %>% filter(width < 6000)
d_neg <- ggviolin(df, x = "exp", y = "width", palette = "Paired", fill = "exp") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    limits = c(0,1.5*(max(df$width)))
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    geom = "pointrange",  # Use pointrange for comparison bars
    position = position_dodge(width = 0.75, height = 10000)  # Adjust the width as needed,
  )+
  # Adding sample size label for each 'exp' group
  geom_text(aes(x = exp, y = 0, label = paste("N =", sample_size)), vjust = -1, hjust = 0.5, size = 4)
ggsave('~/Desktop/test_plot.pdf', d_neg)


write.csv(df,'~/Desktop/test_data.csv')
