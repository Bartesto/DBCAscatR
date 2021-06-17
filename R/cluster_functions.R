# Main clustering and majority functions

#' Ingest of filtered data to list container
#'
#' \code{ingest} reads in the filtered numerical alleles and sample errors,
#'     stores them in a list and converts the same into other formats for
#'     further processing.
#'
#' @details When run it will read in the filtered numerical alleles from the
#'     `results/threshold/` directory and sample errors data from the `results/`
#'     directory. It then stores them in a list that will be accessed by
#'     functions downstream. It also creates a data frame of all possible
#'     combinations of samples and a matrix that will be used for calculating
#'     dissimilarity. These are also stored in the list.
#'
#' @param filtered_alleles File name of the filtered numerical allele data csv
#'     as a character string.
#'
#' @param errors File name of the sample error results data csv as a character
#'     string.
#'
#' @return A list holding 4 data sets. The filtered numerical alleles, the sample
#'     errors, a data frame of sample combinations and a matrix ready for
#'     dissimilarity calculations.
#'
#'@examples
#'\dontrun{
#'data_list <- ingest(filtered_alleles = "numerical_alleles_filtered_at0.8_mt0.2.csv",
#'    errors = "sample_error_results.csv")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import readr
ingest <- function(filtered_alleles, errors){
  suppressWarnings({
    # inputs
    num_out <- readr::read_csv(here::here("results", "threshold", filtered_alleles),
                               col_types = cols())
    results_out <- readr::read_csv(here::here("results", errors),
                                   col_types = cols())

    # names for making matrices
    samp_names <- unique(num_out[['sample']])
    samp_rows <- samp_names[-1]
    samp_cols <- samp_names[-length(samp_names)]

    # output df for unique combos of pair comparisons
    out_df <- data.frame(t(combn(samp_names, 2)), mm = 0, stringsAsFactors = FALSE)

    # output df for conversion to dissimilar matrix
    out_df2 <- data.frame(matrix(NA, nrow = length(samp_names)-1,
                                 ncol = length(samp_names)-1))
    colnames(out_df2) <- samp_cols
    row.names(out_df2) <- samp_rows

    out <- list(out_df, out_df2, num_out, results_out)
    names(out) <- c("combo", "diss", "num_out", "results_out")
    return(out)
  })
}

#' Calculate mismatches
#'
#' \code{calc_mismatches} calculates mismatches between samples and stores
#'     results in the input list.
#'
#' @details When run it calculates the mismatch between all combinations of
#'     samples which have been stored in the input list. It then writes the
#'     results back into the input list.
#'
#' @param slist A list previously compiled by running \code{\link{ingest}}
#'
#' @return An updated list containing calculated mismatches between samples.
#'
#'@examples
#'\dontrun{
#' calc_mismatches(slist = data_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @importFrom  ICSNP pair.diff
calc_mismatches <- function(slist){
  suppressWarnings({
    # inputs
    out_df <- slist[['combo']]
    num_out <- slist[['num_out']]

    # strip off sample
    nums <- num_out[-1]

    # make NAs impossibly big
    nums[is.na(nums)] <- 1000

    # find pairwise diffs
    pdfs <- ICSNP::pair.diff(as.matrix(nums))

    # recode NAs to 0
    pdfs[abs(pdfs) > 500] <- 0

    # recode leftovers to mismatch
    pdfs[pdfs != 0] <- 1

    out_df[['mm']] <- rowSums(pdfs)

    # calculate mismatches between unique pairs - slow
    # allcombo <- length(out_df[,1])
    # for(i in seq_along(out_df[,1])){
    #   cat("missmatch calc ", i , " of ", allcombo, "\n")
    #   samp <- out_df[i, 1]
    #   samp2 <- out_df[i, 2]
    #
    #   df <- num_out %>%
    #     dplyr::filter(sample == samp) %>%
    #     tidyr::pivot_longer(cols = starts_with("x"),
    #                         names_to = "vars",
    #                         values_to = "val")
    #
    #   df2 <- num_out %>%
    #     dplyr::filter(sample == samp2) %>%
    #     tidyr::pivot_longer(cols = starts_with("x"),
    #                         names_to = "vars",
    #                         values_to = "val")
    #
    #   df3 <- df %>%
    #     dplyr::mutate(val2 = df2$val,
    #                   s = val + val2,
    #                   mm = case_when(
    #                     s > 200 & val != val2 ~ 1,
    #                     TRUE ~ 0
    #                   )) %>%
    #     dplyr::summarise(tmm = sum(mm))
    #
    #   out_df[i, 3] <- df3[[1]]
    #
    # }

    # update combos with values
    slist[['combo']] <- out_df
    return(slist)
  })
}

#' Calculate dissimilarity matrix
#'
#' \code{calc_dissimilarity} calculates dissimilarity matrix and stores in the
#'     input list.
#'
#' @details When run it calculates the dissimilarity matrix and writes the
#'     results back into the input list.
#'
#' @param slist1 A list previously compiled by running \code{\link{calc_mismatches}}
#'
#' @return An updated list containing calculated dissimilarity matrix.
#'
#'@examples
#'\dontrun{
#' calc_mismatches(slist = data_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import dplyr
calc_dissimilarity <- function(slist1){
  suppressWarnings({
    # inputs
    out_df <- slist1[['combo']]
    out_df2 <- slist1[['diss']]

    # grab uniques (matches colnames in dissim matrix)
    u_rows <- unique(out_df[,1])

    # populate disimilar df with mismatch values
    for(i in seq_along(u_rows)){
      s <- u_rows[i]

      out <- out_df %>%
        dplyr::filter(X1 == s) %>%
        dplyr::select(mm)

      out_df2[i:length(out_df2[1,]), i] <- out[[1]]
    }

    # update dissimilarity matrix with values
    slist1[['diss']] <- out_df2

    # grab values and names
    tempvect <- as.vector(na.omit(unlist(out_df2)))
    all_names <-  unique(c(colnames(out_df2), row.names(out_df2)))

    # create dist class objects (attributes describe our matrix)
    dist <- structure(tempvect, Size = length(all_names), Labels = all_names,
                      Diag = FALSE, Upper = FALSE, method = "euclidean", #Optional
                      class = "dist")

    distheat <- structure(tempvect, Size = length(all_names), Labels = all_names,
                          Diag = TRUE, Upper = TRUE, method = "euclidean", #Optional
                          class = "dist")

    # add dist class object to main outputs list
    slist1[['dist']] <- dist
    slist1[['distheat']] <- distheat
    return(slist1)
  })
}

#' Calculate mismatches and dissimilarity matrix
#'
#' \code{dissimilarity} calculates mismatches and dissimilarity matrix between
#'     samples from filtered data.
#'
#' @details When run it takes filtered data and calculates mismatches between
#'     samples and a dissimilarity matrix. Results are stored in a list.
#'     Comparing samples for the mismatches takes time and progress status is
#'     printed to screen. Distance objects are also created and can be used to
#'     create visualisations to aid in determining paramaters for group
#'     membership.
#'
#' @inheritParams ingest
#'
#' @return An updated list containing calculated mismatches between samples and
#'     calculated dissimilarity matrix. It will also contain distance objects
#'     that will aid in visualising group memberships.
#'
#'@examples
#'\dontrun{
#'dist <- dissimilarity(filtered_alleles = "numerical_alleles_filtered_at0.8_mt0.2.csv",
#'    errors = "sample_error_results.csv")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import dplyr
#' @import tidyr
#' @import here
#' @import readr
#'
#' @export
dissimilarity <- function(filtered_alleles, errors){
  suppressWarnings({
    slist <- ingest(filtered_alleles, errors)
    slist1 <- calc_mismatches(slist)
    slist2 <- calc_dissimilarity(slist1)
    return(slist2)
  })
}

#' Misassign probability plots and table
#'
#' \code{misassign} assesses different mismatch thresholds by comparing overlap
#'     between allele mismatch distribution of samples assigned to the same group
#'     versus different groups and calculates probability of misassignment of
#'     group membership then summarise them into plots and a table.
#'
#' @details The function compares different mismatch thresholds by generating the
#'     distribution of the pairwise allele mismatch scores for each sample pair.
#'     This distribution is then separated into two groups: allele mismatches
#'     between samples assigned to the same group (i.e. mismatches between scats from
#'     the same putative individual) and allele mismatches between samples assigned
#'     to different groups (i.e. mismatches between scats from different putative
#'     individuals). To assess individual identification success, the same and
#'     different group mismatch distributions are ranked and the upper and lower
#'     0.5 percentiles are calculated. If the difference between the lower and
#'     the upper 0.5 percentile is positive (the overlap column in the summary
#'     table), this means that the distributions are less overlapped and < 1% of
#'     samples have been wrongly assigned. In addition, the probability of
#'     misassignment is  calculated using the `overlap` function in the
#'     package `birdring` with 100,000 simulation and the upper and lower
#'     parameter space set at 99.5% and 0.5% (default).
#'
#'     Outputs of this function generates a series of plots for different thresholds
#'     and a table summary. Each plot consists of the “within” group distribution
#'     in red and “between” groups distribution in blue. The upper 0.5 percentiles
#'     of “within” group distribution and the lower 0.5 percentile of “between”
#'     groups distribution are plotted in dash lines. The number of individuals
#'     indicates the total number of groups identified from each threshold. The
#'     probability of misassignment is calculated with the “overlap” function as
#'     described above. The table summary consists of the following columns: h
#'     indicates the threshold number, ind indicates the total number of groups
#'     identified by each threshold, upper shows the upper 0.5 percentile value
#'     of the “within” group distribution, lower shows the lower 0.5 percentile
#'     value of the “between” groups distribution, overlap is the difference
#'     between the upper and lower 0.5 percentiles columns, and
#'     prob_misassign is the probability of misassignment.
#'
#' @inheritParams dendro_plot
#' @param maxh Integer. It is the maximum "height", or in this case allowable
#'     mismatches, and impacts the length of the x-axis. Default is 10.
#' @param lt Numeric. Lower threshold. Represents the lower tail of the between
#'     individuals distribution. Defaults to 0.005.
#' @param ut Numeric. Upper threshold. Represents the upper tail of the within
#'     individuals distribution. Defaults to 0.995.
#'
#' @return A number of histogram plots, governed by maxh parameter, showing
#'     possible overlap and probability of misassignment are written to jpg file.
#'     A csv summary of key values is also saved.
#'
#'@examples
#'\dontrun{
#'misassign(dist = dissimilarity_list, maxh = 5)
#'}
#'
#' @author Rujiporn Thavornkanlapachai, \email{rujiporn.sun@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @importFrom tibble tibble
#' @importFrom dendextend cutree
#' @importFrom dplyr bind_rows filter
#' @importFrom birdring overlap
#' @import ggplot2
#' @import here
#'
#' @export
misassign <- function(dist, maxh = 10, lt = 0.005, ut = 0.995){
  suppressWarnings({
    distobj <- dist[['dist']]
    clust <- hclust(distobj, method = "average")
    # output tibbles
    group_sum <- tibble::tibble()
    misassign_result <- tibble::tibble()
    # group summaries
    for(i in 1:maxh){
      grp <- dendextend::cutree(clust, k = NULL, h = i)
      sum <- tibble::tibble(sample = names(grp), grp = grp, h = i)
      group_sum <- dplyr::bind_rows(group_sum, sum)
    }

    for(i in 1:maxh){
      # Subset each h
      h_ind <- i
      group_sum_sub <- group_sum %>%
        dplyr::filter(h == h_ind)
      # max group
      g_ind <- max(group_sum_sub[['grp']])

      out_df <- dist[['combo']]

      # Need to generate a new table so each loop doesn't go over the old one and merge each h to the genetic distance data frame
      out_df3 <- out_df # pull out out_df from line 109
      out_df3 <- merge(out_df3, group_sum_sub[, c("sample", "grp")], by.x=c("X1"), by.y=c("sample"))
      out_df3 <- merge(out_df3, group_sum_sub[, c("sample", "grp")], by.x=c("X2"), by.y=c("sample"))

      # Add type of comparison (same versus different groups/individuals)
      out_df3$type <- ifelse(out_df3$grp.x == out_df3$grp.y, "same", "diff")

      # Calculate 99.5% and 0.5% percentile
      lower <- quantile(subset(out_df3,type == 'diff')$mm, lt)
      upper <- quantile(subset(out_df3,type == 'same')$mm, ut)
      overlap <- lower - upper # This number should be positive. If negative means overlap.

      # Calculate probability of misassignment
      pm <- birdring::overlap(subset(out_df3,type == 'diff')$mm,
                              subset(out_df3,type == 'same')$mm,
                              from = lt, to = ut, nsim = 100000,
                              edge.of.parameter.space = F)
      prob_misassign <- format(pm, scientific = F)

      # Create a summary data.frame for each h
      run <- tibble::tibble(h = h_ind, groups = g_ind, upper, lower,
                            overlap, prob_misassign)

      # Append each run to the result data frame
      misassign_result <- dplyr::bind_rows(misassign_result, run)

      # Plot
      pname = here::here("results", "cluster",
                         paste0("h", h_ind, "_genetic_distance.png"))

      p1 <- ggplot(out_df3, aes(x=mm, fill=type)) +
        geom_histogram(data = subset(out_df3, type == 'same'), fill = "red",
                       alpha = 0.2, bins = 40) +
        geom_vline(xintercept = upper, colour = "red", linetype = "longdash") +
        geom_histogram(data = subset(out_df3, type == 'diff'), fill = "blue",
                       alpha = 0.2, bins = 40) +
        geom_vline(xintercept = ifelse(lower == upper, lower +0.15, lower),
                   colour="blue", linetype = "longdash") +
        labs(title = paste(" Number of individuals = ", g_ind, "\n",
                           "Probability of misassignment =",
                           prob_misassign),
             x = "Number of mismatched alleles",
             y = "Number of pairwise comparison") +
        theme_classic()

      ggsave(p1, filename = pname)

    }
    #write out table
    readr::write_csv(misassign_result,
                     here::here("results", "cluster", "misassignment_table.csv"))

  })
}

#' Produces a numerical allele data set at a defined h value
#'
#' \code{group_membership} creates a data frame of numerical alleles based on a
#'     provided h value.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and a maximum number of allowable mismatches,
#'     or h value, to assign group membership to each sample.
#'
#' @inheritParams dendro_plot
#' @param h Integer. Represents determined mismatch value.
#' @param filtered_alleles Character vector. Name of the filtered allele data
#'     set created by running \code{\link{miss_threshold}}.
#'
#' @return A data frame object of numerical alleles with group assignation to be
#'     used in further processing also a csv file ("...withGroups.csv") written
#'     to the `results/cluster/` sub-directory.
#'
#'@examples
#'\dontrun{
#'group_data <- group_membership(dist = dissimilarity_list, h = 5,
#'filtered_alleles = "numerical_alleles_filtered_at0.8_mt0.2.csv")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
#' @import dplyr
#'
#' @export
group_membership <- function(dist, h, filtered_alleles){
  suppressWarnings({
    distobj <- dist[['dist']]
    clust <- hclust(distobj, method = "average")

    # cut based on provide h
    subgrp <- cutree(clust, k = NULL, h = h)

    # get num_out data
    num_out <- readr::read_csv(here::here("results", "threshold",
                                          filtered_alleles), col_types = cols())

    # add to numerical data
    mmout_df <- num_out %>%
      dplyr::mutate(group = subgrp) %>%
      dplyr::select(group, everything()) %>%
      dplyr::arrange(group) %>%
      dplyr::mutate(sample = paste0("ID_", sample))

    # write to file errors per sample
    readr::write_csv(mmout_df,
                     here::here("results", "cluster",
                                paste0("hclust_numerical_mismatch_h",
                                       h, "_withGroups.csv")))
    return(mmout_df)
  })
}

#' Calculates the group majorities and ties data.
#'
#' \code{majorities} creates a csv table of numerical alleles with group
#'     assignation and majority vote values with indications of any within group
#'     ties.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}}, a maximum number of allowable mismatches
#'     (or h value), the filtered allele data set from running
#'     \code{\link{miss_threshold}} and the sample errors from running
#'     \code{\link{gen_errors}} to create a csv. The output contains the numerical
#'     alleles with group assignation in a wide format which also displays the
#'     majority vote per group and any ties.
#'
#' @inheritParams group_membership
#' @param errors Character vector. The name of the sample errors csv.
#'
#' @return A csv of numerical alleles with group assignation, majority vote and
#'     indicated ties written to the `results/cluster/` sub-directory.
#'
#'@examples
#'\dontrun{
#'majorities(dist = dissimilarity_list, h = 5,
#'filtered_alleles = "numerical_alleles_filtered_at0.8_mt0.2.csv",
#'errors = "sample_error_results.csv")
#'}
#'
#'@author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
#' @import dplyr
#' @importFrom tidyr pivot_wider pivot_longer
#' @import here
#'
#' @export
majorities <- function(dist, h, filtered_alleles, errors){
  suppressWarnings({
    mmout_df <- group_membership(dist, h, filtered_alleles)

    results_out <- readr::read_csv(here::here("results", errors),
                                   col_types = cols())

    # function to return mode of a vector with NAs
    Mode <- function(x, na.rm = FALSE) {
      if(na.rm){
        x = x[!is.na(x)]
      }
      ux <- unique(x)
      return(ux[which.max(tabulate(match(x, ux)))])
    }

    # function to indicate "tie" status
    TieMode <- function(x, na.rm = FALSE) {
      if(na.rm){
        x = x[!is.na(x)]
      }
      ux <- unique(x)
      res <- tabulate(match(x, ux))
      tie <- ifelse(sum(res %in% max(res)) > 1, "tie", "ok")
      return(tie)
    }

    # calculate majority vote (using mode)
    d1 <- mmout_df %>%
      tidyr::pivot_longer(cols = starts_with("x"),
                          names_to = "vars",
                          values_to = "val") %>%
      dplyr::group_by(group, vars) %>%
      dplyr::summarise(maj = Mode(val, na.rm = TRUE)) %>%
      tidyr::pivot_wider(names_from = vars,
                         values_from = maj) %>%
      dplyr::mutate(sample = "majority") %>%
      dplyr::select(group, sample, starts_with("x")) %>%
      dplyr::bind_rows(mmout_df) %>%
      dplyr::arrange(group, sample) %>%
      dplyr::ungroup()

    d2 <-  mmout_df %>%
      tidyr::pivot_longer(cols = starts_with("x"),
                          names_to = "vars",
                          values_to = "val") %>%
      dplyr::group_by(group, vars) %>%
      dplyr::summarise(tie = TieMode(val, na.rm = TRUE)) %>%
      dplyr::mutate(tie2 = ifelse(tie == "tie", 1, 0)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-tie) %>%
      tidyr::pivot_wider(names_from = vars,
                         values_from = tie2) %>%
      dplyr::mutate(sample = "tie") %>%
      dplyr::select(group, sample, starts_with("x"))

    d3 <- dplyr::bind_rows(d1, d2) %>%
      dplyr::arrange(group, sample)

    avg_sam_amp_rate <- results_out[,1:2] %>%
      dplyr::mutate(sample = gsub("s", "ID_", sample))

    d4 <- dplyr::left_join(d3, avg_sam_amp_rate, by = "sample")
    #write out table
    readr::write_csv(d4,
                     here::here("results", "cluster",
                                paste0("hclust_numerical_mismatch_h", h,
                                       "_group_majorities_and_ties.csv")))
  })

}

#' Creates a formatted html table of the majorities data export.
#'
#' \code{majorities_html} turns the majorities export csv into a formatted html
#'     table for easy querying.
#'
#' @details Takes the csv output from \code{\link{majorities}} and re-formats it
#'     to an html table. Any ties can be queried thereby guiding manual edits of
#'     the final groups data.
#'
#' @param majorities_csv Character vector. Name of the majorities and ties csv.
#'
#' @return An html table of numerical alleles with group assignation, majority
#'     vote and indicated ties written to the `results/cluster/` sub-directory.
#'
#'@examples
#'\dontrun{
#'majorities_html(majorities_csv = "hclust_numerical_mismatch_h4_group_majorities_and_ties.csv")
#'}
#'
#'@author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import here
#' @import dplyr
#' @importFrom  kableExtra cell_spec kable kable_styling scroll_box save_kable
#' @import readr
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_wider pivot_longer
#'
#' @export
majorities_html <- function(majorities_csv){
  suppressWarnings({
    d <- readr::read_csv(here::here("results", "cluster", majorities_csv),
                         col_types = cols())

    grps <- unique(d[['group']])
    tbl_df <- tibble::tibble()

    # function for bulk of conditional formatting
    diff_to_maj <- function(x){
      ifelse(x == d1$majority,
             kableExtra::cell_spec(x, "html", color = "green", italic = T),
             kableExtra::cell_spec(x, "html", color = "red", bold = T))
    }

    # loop over groups and apply conditional formatting
    for(i in seq_along(grps)){
      d1 <- d %>%
        dplyr::filter(group == grps[i]) %>%
        dplyr::filter(sample != "tie") %>%
        dplyr::select(-avg_amp_rate) %>%
        tidyr::pivot_longer(cols = starts_with("x"), names_to = "vars",
                            values_to = "vals") %>%
        tidyr::pivot_wider(names_from = sample, values_from = vals) %>%
        dplyr::mutate(across(everything(), ~replace_na(.x, 0)))

      d2 <- d1 %>%
        dplyr::mutate_at(vars(starts_with("ID_")), diff_to_maj) %>%
        dplyr::mutate(majority = as.character(majority)) %>%
        dplyr::mutate(majority = kableExtra::cell_spec(majority, background = "darkgrey",
                                                       color = "white", align = "center")) %>%
        tidyr::pivot_longer(cols = c(-group, -vars), names_to = "sample",
                            values_to = "vals") %>%
        tidyr::pivot_wider(names_from = vars, values_from = vals)

      tbl_df <- dplyr::bind_rows(tbl_df, d2)
    }

    # make and view the table. Export as html from the plot window if needed
    tbl_df %>%
      kableExtra::kable(format = "html", escape = F) %>%
      kableExtra::kable_styling("striped", fixed_thead = T, full_width = F) %>%
      kableExtra::scroll_box(width = "100%", height = "800px") %>%
      kableExtra::save_kable(file = here::here("results", "cluster",
                                               "majorities_table.html"),
                             self_contained = TRUE)
  })
}

#' Summary tables with additional site and spatial information
#'
#' \code{summary_tables} joins individuals to site location metadata to create
#'     summary tables.
#'
#' @details Takes the groups csv, from running \code{\link{majorities}} and a
#'     provided lookup table of metadata per sample and provides three summary
#'     tables:
#'     \itemize{
#'       \item capture history
#'       \item individual by date and site
#'       \item site stats
#'     }
#'
#'     The user is required to map required fields from the lookup table to the
#'     function parameters. At a minimum the lookup should contain the
#'     following:
#'     \itemize{
#'       \item sample ID
#'       \item site ID
#'       \item collection date
#'       \item latitude (WGS84)
#'       \item longitude (WGS84)
#'     }
#'
#' @param groups_csv Character vector. Name of the numerical allele groups csv.
#' @param metadata Character vector. Name of the lookup csv which should be
#'     located in the `source/` directory.
#' @param prefix Character vector. Sample names in csv outputs often have a prefix
#'     such as "ID_" to avoid starting with a numeral and the lookup data often
#'     wont. Provide the prefix in the calculated data so that it can be matched
#'     to the lookup data.
#' @param sample Character vector. Name of the column in the lookup data that
#'     contains the sample IDs.
#' @param site_ID Character vector. Name of the column in the lookup data that
#'     contains the site/cave IDs.
#' @param field_date Character vector. Name of the column in the lookup data that
#'     contains the sampling date.
#' @param lat Character vector. Name of the column in the lookup data that
#'     contains the latitude (WGS84) of the sample collection site.
#' @param long Character vector. Name of the column in the lookup data that
#'     contains the longitude (WGS84) of the sample collection site.
#'
#' @return Three csv tables are written to the `results/finalised/` sub-directory
#'     as described in details above.
#'
#'@examples
#'\dontrun{
#'summary_tables(groups_csv = "hclust_numerical_mismatch_h4_withGroups.csv",
#'metadata = "lookup.csv", prefix = "ID_", sample = "sample", site_ID = "roost_name",
#'field_date = "collection_date", lat = "dec_lat", long = "dec_long")
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
#' @importFrom stringr str_split
#' @import lubridate
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble tibble
#' @importFrom zoo as.yearmon
#'
#' @export
summary_tables <- function(groups_csv, metadata, prefix, sample, site_ID, field_date, lat, long){
  suppressWarnings({
    # standardise col names
    rename_cols <- function(site_ID, field_date, lat, long){
      d_clean <- readr::read_csv(here::here("source", metadata),
                                 col_types = cols()) %>%
        dplyr::select(!!sample, !!site_ID, !!field_date, !!lat, !!long)
      names(d_clean) <- c("sample", "site", "date", "lat", "long")
      return(d_clean)
    }
    ngrps <- stringr::str_split(groups_csv, pattern = "_")[[1]][4]
    d1 <- rename_cols(site_ID, field_date, lat, long) %>%
      dplyr::mutate(date = lubridate::parse_date_time(date, c("dmY", "ymd")),
                    date = lubridate::as_date(date),
                    year = lubridate::year(date),
                    month = lubridate::month(date, label = TRUE))

    d2 <- readr::read_csv(here::here("results", "cluster", groups_csv),
                          col_types = cols()) %>%
      dplyr::select(group, sample) %>%
      dplyr::mutate(sample = gsub(pattern = prefix, "", sample)) %>%
      dplyr::left_join(d1, by = "sample")

    ## tables by individuals

    # construct full range of possible dates
    mindate <- min(d2[['date']], na.rm = TRUE)
    maxdate <- max(d2[['date']], na.rm = TRUE)

    dum_dates <- tibble::tibble(date = seq(mindate, maxdate, by = "month")) %>%
      dplyr::mutate(ym = zoo::as.yearmon(date)) %>%
      dplyr::select((-date))

    # first table count of scats by individual, per site, per month
    t1 <- d2 %>%
      dplyr::select(-sample, - lat, -long, -year, -month) %>%
      dplyr::group_by(group, site, date) %>%
      dplyr::count() %>%
      dplyr::rename(scats = "n") %>%
      dplyr::arrange(date) %>%
      dplyr::mutate(date = zoo::as.yearmon(date)) %>%
      tidyr::pivot_wider(names_from = date, values_from = scats) %>%
      tidyr::pivot_longer(cols = -c(group, site), names_to = "ym", values_to = "scats") %>% #, total_scats
      dplyr::filter(!is.na(scats)) %>%
      tidyr::pivot_wider(names_from = site, values_from = scats) %>%
      dplyr::arrange(group) %>%
      dplyr::ungroup() %>%
      dplyr::rowwise(group, ym) %>%
      dplyr::mutate(total_scats = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
      dplyr::rename(date = ym,
                    individual = group,
                    `total scats` = total_scats)

    readr::write_csv(t1,
                     file = here::here("results",
                                       "finalised",
                                       paste0("individual_by_site_by_month_",
                                              ngrps, ".csv")))

    # second table minus length in months
    t2 <- d2 %>%
      dplyr::select(-sample, - lat, -long, -year, -month) %>%
      dplyr::mutate(ym = zoo::as.yearmon(date)) %>%
      dplyr::full_join(dum_dates, by = "ym") %>%
      dplyr::select(-date) %>%
      dplyr::group_by(group, site, ym) %>%
      dplyr::count() %>%
      dplyr::rename(scats = "n") %>%
      dplyr::arrange(ym) %>%
      dplyr::mutate(scats = ifelse(is.na(group), 0, scats)) %>%
      tidyr::pivot_wider(names_from = ym, values_from = scats) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(group)) %>%
      dplyr::rowwise(group, site) %>%
      dplyr::mutate(total_scats = sum(c_across(where(is.numeric)), na.rm = TRUE))

    # calculate total months (add one later)
    t2_1 <- d2 %>%
      dplyr::group_by(group, site) %>%
      dplyr::summarise(minmth = min(date),
                       maxmth = max(date)) %>%
      dplyr:: mutate(int = lubridate::interval(minmth, maxmth),
                     total_mths = round(lubridate::time_length(int, "month"), 0)) %>%
      dplyr::select(group, site, total_mths)

    # second table capture history
    t3 <- t2 %>%
      dplyr::left_join(t2_1, by = c("group", "site")) %>%
      dplyr::mutate(total_mths = ifelse(total_mths != 0, total_mths + 1,
                                        total_mths)) %>%
      dplyr::rename(individual = group,
                    `total scats` = total_scats,
                    `total months` = total_mths) %>%
      dplyr::arrange(individual)
    readr::write_csv(t3, file = here::here("results",
                                           "finalised",
                                           paste0("capture_history_",
                                                  ngrps, ".csv")))


    ## table by site
    t4 <- d2 %>%
      dplyr::select(-sample, - lat, -long, -year, -month) %>%
      distinct() %>%
      dplyr::group_by(site, date) %>%
      dplyr::count() %>%
      dplyr::rename(`# individuals` = "n") %>%
      dplyr::mutate(date = zoo::as.yearmon(date))

    readr::write_csv(t4, file = here::here("results",
                                           "finalised",
                                           paste0("site_history_",
                                                  ngrps, ".csv")))
  })

}

#' Converts majorities output to structure format
#'
#' \code{structure_format} converts majorities output to structure format.
#'
#' @details Takes the majorities and ties csv, from running
#'     \code{\link{majorities}} and converts to a "structure" format which is
#'     saved as a csv to the `results/finalised/` sub-directory.
#'
#' @inheritParams majorities_html
#'
#' @return a csv in "structure" format of the majority data.
#'
#'@examples
#'\dontrun{
#'structure_format(majorities_csv = "hclust_numerical_mismatch_h4_group_majorities_and_ties.csv")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
#' @importFrom stringr str_split
#' @import dplyr
#' @import here
#'
#' @export
structure_format <- function(majorities_csv){
  suppressWarnings({
    d <- readr::read_csv(here::here("results", "cluster", majorities_csv),
                         col_types = cols())

    # grab only majorities, reformat and rearrange
    structure_format <- d %>%
      dplyr::filter(sample == "majority") %>%
      dplyr::mutate(group = paste0("Group", group)) %>%
      dplyr::select(-avg_amp_rate)

    # convert NA to -9
    structure_format[is.na(structure_format)] <- -9

    # get rid of first column name to match structure format
    # doesnt stay like this when written to file
    # colnames(structure_format)[1] <- ""


    oname <- paste0(paste0(stringr::str_split(majorities_csv, "_")[[1]][1:6],
                           collapse = "_"), "_structure_format.csv")
    readr::write_csv(structure_format, file = here::here("results",
                                                         "finalised",
                                                         oname))
  })
}
