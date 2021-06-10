# Main plotting functions

#' Amplification Plots
#'
#' \code{amp_splots} produces amplification plots from the sample error results.
#'
#' @details When run it uses the sample error results, an output from running
#'     \code{\link{gen_errors}}, to produce a panel plot to aid in visualising
#'     amplification rate vs allele error and allelic drop out.
#'
#'     The plots are printed to screen and also saved as a jpg file to the
#'     `results/threshold/` sub-directory. It can be used to determine an amplification
#'     threshold for further data cleaning. Will error with message if no errors
#'     in raw data.
#'
#' @return It will write to jpg file a panel plot of amplification rates.
#'
#' @examples
#' \dontrun{
#' amp_splots()
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom readr read_csv
#' @import here
#' @import ggplot2
#' @importFrom cowplot plot_grid
#'
#' @export
amp_splots <- function(){
  suppressWarnings({
    # grab data
    dat <- readr::read_csv(here::here("results", "sample_error_results.csv"),
                           col_types = cols())

    # check dimensions logic no duplicate (no errors) data
    if(dim(dat)[2] > 2){
      error_out <- dat %>%
        tidyr::drop_na(allele_error)

      # plot amplification rate vs allele error
      p1 <- ggplot(error_out, aes(x = avg_amp_rate, y = allele_error)) +
        geom_point(size = 2) +
        geom_smooth(method = lm, se = FALSE) +
        labs(title = "Amplification rate vs allele error",
             x = "average amplification rate",
             y = "allele error") +
        theme_bw()

      # plot amplification rate vs allelic dropout
      p2 <- ggplot(error_out, aes(x = avg_amp_rate, y = allelic_drop_out)) +
        geom_point(size = 2) +
        geom_smooth(method = lm, se = FALSE) +
        labs(title = "Amplification rate vs allelic drop out",
             x = "average amplification rate",
             y = "allelic drop out") +
        theme_bw()

      p3 <- cowplot::plot_grid(p1, p2, labels = "auto", ncol = 1)

      fname <- here::here("results", "threshold", "amplification_v_allele_scatterplots.jpg")
      ggsave2(fname, p3)
      print(p3)
    } else {
      stop(paste0("Not multivariate. Data probably has no error values and was ",
                  "most likely generated from raw data with no replicates"))
    }
  })

}

#' Data "missingness" plot
#'
#' \code{miss_hist} creates a histogram to visualise NA proportions for the loci.
#'
#' @details When run it uses the "loci_NA.csv", an output
#'     from running \code{\link{amp_threshold}}, to produce a histogram to aid in
#'     visualising data "missingness".
#'
#'     The plot is printed to screen and also saved as a jpg file to
#'     the `results/threshold/` sub-directory. It can be used to determine a
#'     "missingness" threshold for further data cleaning.
#'
#' @return It will write to jpg file a histogram of NA proportions for the loci.
#'
#' @examples
#' \dontrun{
#' miss_hist()
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom readr read_csv
#' @import here
#' @import ggplot2
#'
#' @export
miss_hist <- function(){
  suppressWarnings({
    dat <- readr::read_csv(here::here("results", "threshold", "loci_NA.csv"),
                           col_types = cols()) %>%
      dplyr::mutate(id = ifelse(stringr::str_detect(loci, "_a"), "a", "b")) %>%
      dplyr::filter(id == "a")

    p1 <- ggplot(dat) +
      geom_histogram(aes(proportion), binwidth = 0.005) +
      labs(title = "Missingness histogram",
           x = "proportion of NAs",
           y = "count") +
      theme_bw()
    ggsave(here::here("results", "threshold", "missingness_histogram.jpg"), p1)
    print(p1)
  })
}

#' Amplification histogram
#'
#' \code{amp_hist} creates a histogram to visualise distribution of average
#'     amplification rates.
#'
#' @details When run it uses the "sample_error_results.csv", an output
#'     from running \code{\link{gen_errors}}, to produce a histogram to aid in
#'     visualising the distribution of amplification.
#'
#'     The plot is printed to screen and also saved as a jpg file to
#'     the `results/threshold/` sub-directory. It can be used to determine an
#'     amplification threshold for further data cleaning.
#'
#' @return It will write to jpg file a histogram average amplification rates.
#'
#' @examples
#' \dontrun{
#' amp_hist()
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom readr read_csv
#' @import here
#' @import ggplot2
#'
#' @export
amp_hist <- function(){
  suppressWarnings({
    dat <- readr::read_csv(here::here("results", "sample_error_results.csv"),
                           col_types = cols()) %>%
      dplyr::mutate(amp = avg_amp_rate)

    p1 <- ggplot(dat) +
      geom_histogram(aes(amp)) +
      labs(title = "Amplification histogram",
           x = "average amplification rate",
           y = "count") +
      theme_bw()
    ggsave(here::here("results", "amplification_histogram.jpg"), p1)
    print(p1)
  })
}

#' Dendrogram plot
#'
#' \code{dendro_plot} creates and dendrogram plot.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates a dendrogram. A dendrogram can
#'     aid in visualising a "height" parameter and corresponding group
#'     membership of samples. Here the "height" parameter equates to number of
#'     mismatches. The plot is written to the `results/cluster/`
#'     sub-directory.
#'
#' @param dist A list object created after running \code{\link{dissimilarity}}.
#'
#' @return It will write to jpg file a dendrogram plot.
#'
#' @examples
#' \dontrun{
#' dendro_plot(dist = dissimilarity_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @import ggplot2
#' @importFrom ggdendro ggdendrogram
#' @import here
#'
#' @export
dendro_plot <- function(dist){
  suppressWarnings({
    distobj <- dist[['dist']]
    clust <- hclust(distobj, method = "average")
    p1 <- ggdendro::ggdendrogram(clust, rotate = TRUE)
    ggsave(here::here("results", "cluster", "dendrogram.jpg"), p1)
    print(p1)
  })
}

#' Elbow plot
#'
#' \code{elbow_plot} creates an elbow plot.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates an elbow plot. An elbow plot can
#'     aid in visualising different "height" parameters and corresponding number
#'     of groups created.  Here the "height" parameter equates to number of
#'     mismatches. The plot is written to the `results/cluster/`
#'     sub-directory.
#'
#' @inheritParams dendro_plot
#' @param maxh is the maximum "height", or in this case allowable mismatches,
#'     as an integer and impacts the length of the x-axis. Default is 10.
#'
#' @return It will write to jpg file an elbow plot.
#'
#' @examples
#' \dontrun{
#' elbow_plot(dist = dissimilarity_list, maxh = 5)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @import ggplot2
#' @import here
#'
#' @export
elbow_plot <- function(dist, maxh = 10){
  suppressWarnings({
    distobj <- dist[['dist']]
    clust <- hclust(distobj, method = "average")
    plotdf <- tibble::tibble()
    for(i in 0:maxh){
      sub_grp <- cutree(clust, k = NULL, h = i)
      df <- data.frame(threshold = i, groups = length(table(sub_grp)))
      plotdf <- dplyr::bind_rows(plotdf, df)
    }
    p1 <- ggplot(plotdf) +
      geom_line(aes(x = threshold, y = groups)) +
      scale_x_continuous(breaks = c(0:maxh)) +
      labs(y = "Groups",
           x = "Mismatch Threshold") +
      theme_bw()
    ggsave(here::here("results", "cluster", "elbow.jpg"), p1)
    print(p1)
  })

}

#' Mismatch histogram
#'
#' \code{freq_hist} creates a histogram of the frequency of mismatches.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates a histogram. The histogram can
#'     help visualise the frequency of mismatches and can aid in selecting a
#'     suitable value for creating groups. The plot is written to the
#'     `results/cluster/` sub-directory.
#'
#' @inheritParams dendro_plot
#'
#' @return It will write to png file a histogram of mismatch frequency.
#'
#' @examples
#' \dontrun{
#' freq_hist(dist = dissimilarity_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @import here
#'
#' @export
freq_hist <- function(dist){
  suppressWarnings({
    #func for max # missmatch
    my_max <- function(x) ifelse(!all(is.na(x)), max(x, na.rm = TRUE), NA)
    dissobj <- dist[['diss']]
    freq <- unlist(dissobj)
    hist_out <- here::here("results", "cluster", "hclust_group_mismatch_histogram.png")
    png(hist_out)
    hist(freq,
         xlab = "No. of mismatch",
         ylab = "Frequency",
         labels = TRUE,
         breaks = -1:my_max(freq))
    dev.off()
  })
}

#' Mismatch heatplot
#'
#' \code{heat_plot} creates a heat plot based on mismatches.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates a heat map. The heat map can
#'     help visualise groupings of samples based on mismatch values. The plot is
#'     written to the `results/cluster/` sub-directory. It also generates a
#'     "plotly" visualisation in the viewer if using RStudio whereby values can
#'     be interrogated by hovering the mouse cursor.
#'
#' @inheritParams dendro_plot
#'
#' @return It will write to png file a heat map.
#'
#' @examples
#' \dontrun{
#' heat_plot(dist = dissimilarity_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom heatmaply heatmaply
#' @importFrom webshot install_phantomjs
#' @import here
#'
#' @export
heat_plot <- function(dist){
  suppressMessages({
    webshot::install_phantomjs()
    distobj <- dist[['dist']]
    m <- as.matrix(distobj)
    heatmaply::heatmaply(m, k_row = 2, k_col = 2,
                         file = here::here("results", "cluster", "hclust_heatmap_plot.html"),
                         width = 1200, height = 750, plot_method = "ggplot")
    heatmaply::heatmaply(m, k_row = 2, k_col = 2,
                         file = here::here("results", "cluster", "hclust_heatmap_plot.png"),
                         width = 1200, height = 750, plot_method = "ggplot")
  })
}

#' Sample correlation heatplot
#'
#' \code{heat_cor_plot} creates a heat plot based on correlations.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates a correlation heat map. The heat
#'     map can help visualise groupings of samples based on correlations of
#'     mismatch values. The plot is written to the `results/cluster/`
#'     sub-directory. It also generates a "plotly" visualisation in the viewer
#'     if using RStudio whereby values can be interrogated by hovering the mouse
#'     cursor.
#'
#' @inheritParams dendro_plot
#'
#' @return It will write to png file a correlation heat map.
#'
#' @examples
#' \dontrun{
#' heat_cor_plot(dist = dissimilarity_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom heatmaply heatmaply_cor
#'
#' @import here
#' @importFrom webshot install_phantomjs
#' @export
heat_cor_plot <- function(dist){
  suppressMessages({
    webshot::install_phantomjs()
    distobj <- dist[['dist']]
    m <- as.matrix(distobj)
    heatmaply::heatmaply_cor(m, k_row = 2, k_col = 2,
                             file = here::here("results", "cluster",
                                               "hclust_heatmap_cor_plot.html"),
                             width = 1200, height = 750, plot_method = "ggplot")
    heatmaply::heatmaply_cor(m, k_row = 2, k_col = 2,
                             file = here::here("results", "cluster",
                                               "hclust_heatmap_cor_plot.png"),
                             width = 1200, height = 750, plot_method = "ggplot")
  })
}

#'
