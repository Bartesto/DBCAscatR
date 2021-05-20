# Main DBCAscatR plotting functions

#' Amplification Plots
#'
#' \code{amp_splots} produces amplification plots from the sample error results.
#'
#' @details When run it uses the sample error results, an output from running
#'     \code{\link{gen_errors}}, to produce a panel plot to aid in visualising
#'     amplification rate vs allele error and allelic drop out.
#'
#'     The plots are printed to screen and also saved as a jpg file to the
#'     `results/` sub-directory. It can be used to determine an amplification
#'     threshold for further data cleaning. Will error with message if no errors
#'     in raw data.
#'
#' @return It will write to jpg file a panel plot of amplification rates.
#'
#' @examples
#' \dontrun{
#' amp_plots()
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

    fname <- here::here("results", "amplification_v_allele_scatterplots.jpg")
    ggsave2(fname, p3)
    print(p3)
  } else {
    stop(paste0("Not multivariate. Data probably has no error values and was ",
                "most likely generated from raw data with no replicates"))
  }

}

#' Data missingness plot
#'
#' \code{miss_plot} creates a histogram to visualise NA proportions for the loci.
#'
#' @details When run it uses the "sample error results "loci_NA.csv", an output
#'     from running \code{\link{amp_threshold}}, to produce a histogram to aid in
#'     visualising data "missingness".
#'
#'     The plot is printed to screen and also saved as a jpg file to the
#'     `results/` sub-directory. It can be used to determine an missingness
#'     threshold for further data cleaning.
#'
#' @return It will write to jpg file a histogram of NA proportions for the loci.
#'
#' @examples
#' \dontrun{
#' miss_plot()
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
miss_plot <- function(){
  suppressWarnings({
    dat <- readr::read_csv(here::here("results", "loci_NA.csv"),
                           col_types = cols()) %>%
      dplyr::mutate(id = ifelse(stringr::str_detect(loci, "_a"), "a", "b")) %>%
      dplyr::filter(id == "a")

    p1 <- ggplot(dat) +
      geom_histogram(aes(proportion), binwidth = 0.005) +
      labs(title = "Missingness histogram",
           x = "proportion of NAs",
           y = "count") +
      theme_bw()
    ggsave(here::here("results", "missingness_histogram.jpg"), p1)
    print(p1)
  })
}
