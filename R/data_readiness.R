# Main data cleaning functions

#' Creates a directory structure for inputs and outputs for a `ScatMatch` workflow.
#'
#' \code{workspace} creates a directory structure for use with the `ScatMatch`
#'     functions.
#'
#' @details When run within an `RStudio` project folder structure (highly
#'     recommended), it will create  `results/` and `source/` sub-directories.
#'     After running, the USER should place a copy of the raw data input csv file to the
#'     newly created `source/` sub-directory. All outputs from `ScatMatch` functions
#'     will be written to the `results/` sub-directory.
#'
#' @return Two sub-directories, `results/` and `source/` will be created in the
#'     USERs `RStudio` project directory.
#'
#'@examples
#'\dontrun{
#'workspace()
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#'
#' @export
workspace <- function(){
  s <- here::here("source")
  r <- here::here("results")
  t_out <- here::here("results", "threshold")
  c_out <- here::here("results", "cluster")
  f_out <- here::here("results", "finalised")
  ifelse(!dir.exists(s),
         dir.create(s), FALSE)
  ifelse(!dir.exists(r),
         dir.create(r), FALSE)
  ifelse(!dir.exists(t_out),
         dir.create(t_out), FALSE)
  ifelse(!dir.exists(c_out),
         dir.create(c_out), FALSE)
  ifelse(!dir.exists(f_out),
         dir.create(f_out), FALSE)
}

#' Ingest of raw data and creation of two tidy data sets.
#'
#' \code{data_in} reads in raw data SNP genotype results file and creates two
#'     clean data sets, one with explicit NA's and the other with replacement
#'     blank cells.
#'
#' @details When run, it reads in a raw data csv file from the `source/` sub-directory
#'     performs some data cleaning and standardisation, and stores two data frames
#'     in a list. One data frame contains explicit NA's whilst the other has NA's
#'     coded as blank cells.
#'
#'     The raw data expected is SNP genotypes generated through the Agena
#'     Bioscience MassARRAY genotyping system as provided by the Australian
#'     Genome Research Facility (AGRF).
#'
#' @param filename File name of raw data input csv as a character string.
#'
#' @param replicates Logical. TRUE if replicates are present in the raw data.
#'
#' @param suffix File suffix to denote sample replicates. Default is the character
#'     string"_dup", change if required. Ignore if replicates = FALSE.
#'
#' @return When assigned to an object it will create a list containing two data
#'     frames.
#'
#'@examples
#'\dontrun{
#'data_in(filename = "CAGRF20021407_raw.csv", suffix = "_dup")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import readr
#' @importFrom janitor clean_names
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
data_in <- function(filename, replicates = TRUE, suffix = "_dup"){
  suppressWarnings({
    if(replicates == TRUE){
      file <- here::here("source", filename)
      dat <- readr::read_csv(file, col_types = cols())
      # logic to catch misnamed replicate
      if(any(stringr::str_detect(dat[['X1']], suffix))){
        # data with NA
        file <- here::here("source", filename)
        na_dat <- readr::read_csv(file, na = "Fail", col_types = cols()) %>%
          janitor::clean_names() %>%
          dplyr::rename(sample = .data$x1) %>%
          dplyr::mutate(rep = case_when(
            stringr::str_detect(sample, suffix) ~ 2,
            TRUE ~ 1),
            sample = gsub(suffix, "", sample),
            ind = stringr::str_detect(sample, "Blank")) %>%
          dplyr::filter(.data$ind != TRUE) %>%
          dplyr::select(sample, rep, everything(), -.data$ind)

        # data with blanks for NA
        all_dat <- read_csv(file, col_types = cols()) %>%
          janitor::clean_names() %>%
          dplyr::rename(sample = .data$x1) %>%
          dplyr::mutate(rep = case_when(
            stringr::str_detect(sample, suffix) ~ 2,
            TRUE ~ 1),
            sample = gsub(suffix, "", sample),
            ind = stringr::str_detect(sample, "Blank")) %>%
          dplyr::filter(.data$ind != TRUE) %>%
          dplyr::select(sample, rep, everything(), -.data$ind) %>%
          tidyr::pivot_longer(cols = starts_with("x"),
                              names_to = "marker",
                              values_to = "val") %>%
          dplyr::mutate(val = ifelse(.data$val == "Fail", "", .data$val)) %>%
          tidyr::pivot_wider(names_from = .data$marker,
                             values_from = .data$val)

        out <- list(na_dat = na_dat, all_dat = all_dat)
        return(out)
      } else {
        stop("The suffix chosen is not present in the raw data")
      }
    } else {
      # data with NA
      file <- here::here("source", filename)
      na_dat <- readr::read_csv(file, na = "Fail", col_types = cols()) %>%
        janitor::clean_names() %>%
        dplyr::rename(sample = .data$x1) %>%
        dplyr::mutate(rep = case_when(
          stringr::str_detect(sample, suffix) ~ 2,
          TRUE ~ 1),
          sample = gsub(suffix, "", sample),
          ind = stringr::str_detect(sample, "Blank")) %>%
        dplyr::filter(.data$ind != TRUE) %>%
        dplyr::select(sample, rep, everything(), -.data$ind)

      # data with blanks for NA
      all_dat <- read_csv(file, col_types = cols()) %>%
        janitor::clean_names() %>%
        dplyr::rename(sample = .data$x1) %>%
        dplyr::mutate(rep = case_when(
          stringr::str_detect(sample, suffix) ~ 2,
          TRUE ~ 1),
          sample = gsub(suffix, "", sample),
          ind = stringr::str_detect(sample, "Blank")) %>%
        dplyr::filter(.data$ind != TRUE) %>%
        dplyr::select(sample, rep, everything(), -.data$ind) %>%
        tidyr::pivot_longer(cols = starts_with("x"),
                            names_to = "marker",
                            values_to = "val") %>%
        dplyr::mutate(val = ifelse(.data$val == "Fail", "", .data$val)) %>%
        tidyr::pivot_wider(names_from = .data$marker,
                           values_from = .data$val)

      out <- list(na_dat = na_dat, all_dat = all_dat)
      return(out)
    }

  })
}

#' Calculates genotyping error rates
#'
#' \code{main_errors} takes the list generated by \code{\link{data_in}} and
#'     generates genotyping error outputs which are written to csv file.
#'
#' @details When run it calculates genotyping error rates (sample amplification
#'     rate and, when replicate samples provided, allelic dropout and false allele
#'     rates) on the complete raw data set which it receives in the form of a
#'     list from the \code{\link{data_in}} function. Outputs are written to csv
#'     file to the `results/` sub-directory as:
#'     \itemize{
#'         \item numerical_alleles.csv
#'         \item sample_error_results.csv
#'         \item summary_error_results.csv}
#'
#'      These outputs are used in downstream processes and are further refined.
#'
#' @param dl data list. This is a list containing two data frames as output from
#'     the \code{\link{data_in}} function.
#'
#' @return It will write to csv file genetic error file csv's for further
#'     processing and/or evaluation.
#'
#'@examples
#'\dontrun{
#'main_errors(dl)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom stats sd
#' @importFrom rlang .data
#' @import readr
#' @import tidyr
main_errors <- function(dl){
  suppressWarnings({
    na_dat <- dl[['na_dat']]
    all_dat <- dl[['all_dat']]
    # just variables
    df_loci <- na_dat  %>%
      dplyr::select(starts_with("X")) # keep loci columns for calcs

    # clean with blanks
    cl_dat_bl <- all_dat[rowSums(is.na(df_loci)) != ncol(df_loci),]

    # clean with na's
    cl_dat_na <- na_dat[rowSums(is.na(df_loci)) != ncol(df_loci),]

    # anonymous function for merge alleles for grouping step
    f1 <- function(x){trimws(paste(x, collapse = ''))}

    # set up for loop variables stuff to use in loop
    lgth_vars <- dim(df_loci)[2]
    u_samp <- unique(unlist(cl_dat_bl[, 1]))
    samp_names <- paste0("s", u_samp)

    # set up data frames to hold results
    results_out <- tibble::tibble()
    num_out <- tibble::tibble()

    # loop and calculate interim and end results
    for(i in seq_along(u_samp)){
      # using na data
      avg_amp_rate <- cl_dat_na %>%
        dplyr::filter(sample == u_samp[i]) %>%
        tidyr::pivot_longer(cols = starts_with("x"),
                            names_to = "vars",
                            values_to = "vals") %>%
        dplyr::group_by(rep) %>%
        dplyr::summarise(amp_rate = sum(!is.na(.data$vals))/lgth_vars) %>%
        dplyr::summarise(avg_amp_rate = mean(.data$amp_rate))

      # interim data (creating new variables to calculate errors)
      d <- cl_dat_bl %>%
        dplyr::filter(sample == u_samp[i]) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise_all(.funs = f1) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_longer(cols = starts_with("x"),
                            names_to = "vars",
                            values_to = "vals") %>%
        dplyr::mutate(
          f = case_when(
            nchar(.data$vals) == 4 ~ .data$vals,
            TRUE ~ ""),
          c = nchar(.data$f),
          loc_err = ifelse(c == 4 & substr(.data$f, 1, 2) != substr(.data$f, 3, 4),
                           0.5, 0),
          p1 = substr(.data$f, 1, 1),
          p2 = substr(.data$f, 2 ,2),
          p3 = substr(.data$f, 3, 3),
          p4 = substr(.data$f, 4 ,4),
          a_mm1 = ifelse(c == 4 & .data$p1 != .data$p3, 1, 0),
          a_mm2 = ifelse(c == 4 & .data$p2 != .data$p4, 1, 0),
          a_err = (.data$a_mm1 + .data$a_mm2)/4,
          het_only = case_when(
            .data$p1 == .data$p2 & .data$p3 == .data$p4 ~ 0,
            TRUE ~ 1),
          a_drop = case_when(
            .data$het_only == 1 & substr(.data$f, 1, 2) != substr(.data$f, 3, 4) ~ 0.5,
            TRUE ~ 0),
          fa = case_when(
            c == 4 & .data$het_only == 0 & substr(.data$f, 1 ,2) != substr(.data$f, 3 ,4) ~ 0.5,
            TRUE ~ 0))

      # take interim data and use to convert to numerical values for export
      num_vals <- d %>%
        dplyr::mutate(new_var = case_when(
          nchar(.data$vals) == 2 ~ substr(.data$vals, 1, 2),
          nchar(.data$vals) == 4 & substr(.data$vals, 1, 2) == substr(.data$vals, 3, 4) ~ substr(.data$vals, 1, 2),
          TRUE ~ "NA"
        )) %>%
        dplyr::mutate(to_recode = case_when(
          new_var == "NA" ~ 0,
          TRUE ~ 1
        ),
        a1 = ifelse(.data$to_recode == 1, substr(.data$vals, 1, 1), NA),
        b1 = ifelse(.data$to_recode == 1, substr(.data$vals, 2, 2), NA),
        a = case_when(
          a1 == "A" ~ 110,
          a1 == "T" ~ 120,
          a1 == "G" ~ 130,
          a1 == "C" ~ 140
        ),
        b = case_when(
          b1 == "A" ~ 110,
          b1 == "T" ~ 120,
          b1 == "G" ~ 130,
          b1 == "C" ~ 140)
        ) %>%
        # dplyr::select(sample, vars, a, b) %>%
        # tidyr::pivot_longer(cols = c("a", "b"),
        #                     names_to = "locus",
        #                     values_to = "value") %>%
        # dplyr::arrange(vars, locus) %>%
        # tidyr::pivot_wider(names_from = c(vars, locus), values_from = value)
        dplyr::select(.data$sample, .data$vars, .data$a, .data$b) %>%
        tidyr::pivot_longer(cols = c("a", "b"),
                            names_to = "rep",#
                            values_to = "value") %>%
        dplyr::arrange(vars, .data$value) %>%
        dplyr::mutate(locus = rep(c("a", "b"), dim(d)[1])) %>%
        dplyr::select(-rep) %>%
        tidyr::pivot_wider(names_from = c(vars, .data$locus), values_from = .data$value)

      # take interim data and calculate errors for export
      results <- d %>%
        dplyr::summarise(locus_error = sum(.data$loc_err)/sum(c == 4),
                         allele_error = sum(.data$a_err)/sum(c == 4),
                         allelic_drop_out = sum(.data$a_drop)/sum(.data$het_only),
                         false_allele = sum(.data$fa)/(sum(c)/4)) %>%
        dplyr::bind_cols(avg_amp_rate) %>%
        dplyr::mutate(sample = samp_names[i]) %>%
        dplyr::select(.data$sample, .data$avg_amp_rate, .data$allele_error,
                      .data$locus_error, .data$allelic_drop_out,
                      .data$false_allele)

      # bind exports to appropriate result data frames
      results_out <- dplyr::bind_rows(results_out, results)
      num_out <- dplyr::bind_rows(num_out, num_vals)
    }

    # clean out errors if a no duplicate data set
    results_out_clean <- results_out %>%
      tidyr::pivot_longer(cols = -sample,
                          names_to = "error",
                          values_to = "value") %>%
      filter(!is.na(.data$value)) %>%
      tidyr::pivot_wider(names_from = .data$error,
                         values_from = .data$value)

    # take the error results and further summarise
    summaries <- results_out_clean %>%
      tidyr::pivot_longer(cols = -sample,
                          names_to = "error",
                          values_to = "value") %>%
      dplyr::group_by(.data$error) %>%
      dplyr::summarise(avg = mean(.data$value, na.rm = TRUE),
                       se = sd(.data$value, na.rm = TRUE)/sqrt(n()))

    # write to file errors per sample
    readr::write_csv(results_out_clean, here::here("results", "sample_error_results.csv"))

    # write to file summary errors for whole of data run
    readr::write_csv(summaries, here::here("results", "summary_error_results.csv"))

    # write to file numerical version
    readr::write_csv(num_out, here::here("results", "numerical_alleles.csv"))
  })
}

#' Main function wrapper to calculate genotyping error rates
#'
#' \code{gen_errors} takes a raw data file and generates genotyping error outputs
#'     which are written to csv file.
#'
#' @details When run it calculates genotyping errors on the complete raw data set
#'     which it reads in from the `source/` sub-directory. The raw data expected
#'     is SNP genotypes generated through the Agena Bioscience MassARRAY
#'     genotyping system as provided by the Australian Genome Research Facility
#'     (AGRF).
#'
#'     Outputs are written to csv file to the `results/` sub-directory as:
#'     \itemize{
#'         \item numerical_alleles.csv
#'         \item sample_error_results.csv
#'         \item summary_error_results.csv}
#'
#'      The genotyping error rates (sample amplification rate and, when replicate
#'      samples provided, allelic dropout and false allele rates) are used in
#'      downstream processes and are further refined.
#'
#' @inheritParams data_in
#'
#' @return It will write to csv file genotype error file csv's for further
#'     processing and/or evaluation.
#'
#'@examples
#'\dontrun{
#'gen_errors(filename = "CAGRF20021407_raw.csv", suffix = "_dup")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import readr
#' @importFrom janitor clean_names
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @import dplyr
#' @import tidyr
#'
#' @export
gen_errors <- function(filename, replicates = TRUE, suffix = "_dup"){
  suppressWarnings({
    # ingest data and make data list
    dl <- data_in(filename, replicates, suffix)
    # munge and create error list
    main_errors(dl)
  })
}

#' Applies quality filtering of samples by choosing an amplification threshold
#' below which samples are discarded and produces interim filtered results.
#'
#' \code{amp_threshold} takes a threshold value and applies it to the numerical
#'     allele data. Also produces a record of which samples were filtered out and
#'     proportions of NA's per loci.
#'
#' @details The predetermined amplification threshold, ascertained through
#'     visualisation \code{\link{amp_splots}}, is used to filter the numerical
#'     alleles data set and written to csv. The names of all the samples that are
#'     rejected at this threshold are also written to csv file.
#'
#'     Lastly the proportion of missing data (NA’s) per locus is written to csv.
#'     All outputs are written to `results/`.
#'
#' @param at Numeric value of the average amplification rate to filter the
#'     numerical alleles. Values greater than or equal to are retained.
#'
#' @return It will write to `results/` three csv files:
#'     \itemize{
#'         \item filtered alleles data with threshold indicated in the name
#'         \item sample that were filtered at with threshold indicated in the name
#'         \item proportion of NA's for the loci}
#'
#'     It will also print to screen the number of samples filtered out at this
#'     threshold.
#'
#'@examples
#'\dontrun{
#'amp_threshold(at = 0.8)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import readr
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_replace
#' @importFrom rlang .data
#'
#' @export
amp_threshold <- function(at){
  suppressWarnings({
    results_out <- readr::read_csv(here::here("results",
                                              "sample_error_results.csv"),
                                   col_types = cols())
    num_out <- readr::read_csv(here::here("results",
                                          "numerical_alleles.csv"),
                               col_types = cols())
    # filter results list by amplification success rate
    filt_out <- results_out %>%
      dplyr::filter(.data$avg_amp_rate >= at)

    filt_names <- filt_out['sample']
    filt_names[['sample']] <- stringr::str_replace(filt_names[['sample']], "s",
                                                   "")
    # filter num_out dataframe by sample name values in filt_out
    num_out_filt_samps <- num_out[num_out[['sample']] %in% filt_names[['sample']], ]

    # write to file interim filtered alleles on average amplification rate (at)
    cname <- paste0("numerical_alleles_filtered_at", at, ".csv")
    readr::write_csv(num_out_filt_samps, here::here("results", "threshold", cname))

    # for records, get a list of samples that were filtered
    samples_x <- results_out %>%
      dplyr::filter(.data$avg_amp_rate < at) %>%
      dplyr::select(.data$sample, .data$avg_amp_rate)

    #print to screen number of removals
    ls <- length(samples_x[['sample']])
    cat("At an amplification threshold  of", at, ": filtering out", ls, "samples")

    # write to file samples filtered out based on average amplification rate (at)
    c2name <- paste0("samples_filtered_at", at, ".csv")
    readr::write_csv(samples_x, here::here("results", "threshold", c2name))

    # summarise proportion NAs per locus in filtered num_out dataframe
    loc_NA <- num_out_filt_samps %>% summarise_all(~ mean(is.na(.))) %>%
      dplyr::select(-sample) %>%
      tidyr::pivot_longer(everything(), names_to = "loci",
                          values_to = "proportion")

    # write to file props of NA per locus for histogram plotting
    readr::write_csv(loc_NA, here::here("results", "threshold", "loci_NA.csv"))

  })
}

#' Applies quality filtering of loci by applying a maximum amount of missing data
#' allowed per locus above which loci are discarded
#'
#' \code{miss_threshold} takes a "missingness" threshold and applies it to a
#'     filtered numerical allele data set and then writes the results to csv.
#'
#' @details The predetermined "missingness" threshold, ascertained through
#'     visualisation \code{\link{miss_hist}}, is used to further filter the
#'     numerical alleles data set and written to csv to the `results/`
#'     sub-directory.
#'
#'     The amplification threshold parameter is used to select the correct
#'     data set from the `results/` sub-directory that had been created by the
#'     \code{\link{amp_threshold}} function.
#'
#' @param mt Numeric value of the "missingness" threshold to further filter the
#'     numerical alleles. Values less than the threshold are retained.
#'
#' @param at Numeric value of the average amplification threshold used on the
#'     input data. NOTE used to select a pre-made data.
#'
#' @return Writes a numerical alleles data csv, filtered to the `at` and `mt`
#'     thresholds, to the `results/` sub-directory. It will also print to screen
#'     the number of loci that will be filtered out at this "missingness" threshold.
#'
#'@examples
#'\dontrun{
#'miss_threshold(mt = 0.2, at = 0.8)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import readr
#' @import dplyr
#' @importFrom rlang .data
#'
#' @export
miss_threshold <- function(mt, at){
  suppressWarnings({

    # read in locus NA proportion data
    loc_NA_names <- readr::read_csv(here::here("results", "threshold", "loci_NA.csv"),
                                    col_types = cols()) %>%
      dplyr::filter(.data$proportion >= mt) %>%
      dplyr::pull(.data$loci)

    #print to screen number of removals
    ln <- length(loc_NA_names)
    cat("At an missingness threshold of", mt, ": filtering out", ln, "loci")

    # make name of correct amplification filtered data and read in
    num_name <- paste0("numerical_alleles_filtered_at", at, ".csv")
    num_out_filt_samps <-  readr::read_csv(here::here("results", "threshold",
                                                      num_name),
                                           col_types = cols())

    # filter loci with high missingness from num_out dataframe
    num_out_filt_loci <- num_out_filt_samps %>%
      dplyr::select(-all_of(loc_NA_names))

    # write to file numerical version that has been quality filtered for samples
    # and loci
    cname <- paste0("numerical_alleles_filtered_at", at, "_mt", mt, ".csv")
    readr::write_csv(num_out_filt_loci, here::here("results", "threshold", cname))

  })
}
