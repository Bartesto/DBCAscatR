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
#' @examples
#' \dontrun{
#' data_list <- ingest(filtered_alleles = "numerical_alleles_filtered_a0.8_m0.2.csv",
#'     errors = "sample_error_results.csv")
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @import here
#' @importFrom readr read_csv
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
#'  @examples
#' \dontrun{
#' calc_mismatches(slist = data_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @import dplyr
#' @import tidyr
calc_mismatches <- function(slist){
  suppressWarnings({
    # inputs
    out_df <- slist[['combo']]
    num_out <- slist[['num_out']]

    # calculate mismatches between unique pairs
    allcombo <- length(out_df[,1])
    for(i in seq_along(out_df[,1])){
      cat("missmatch calc ", i , " of ", allcombo, "\n")
      samp <- out_df[i, 1]
      samp2 <- out_df[i, 2]

      df <- num_out %>%
        dplyr::filter(sample == samp) %>%
        tidyr::pivot_longer(cols = starts_with("x"),
                            names_to = "vars",
                            values_to = "val")

      df2 <- num_out %>%
        dplyr::filter(sample == samp2) %>%
        tidyr::pivot_longer(cols = starts_with("x"),
                            names_to = "vars",
                            values_to = "val")

      df3 <- df %>%
        dplyr::mutate(val2 = df2$val,
                      s = val + val2,
                      mm = case_when(
                        s > 200 & val != val2 ~ 1,
                        TRUE ~ 0
                      )) %>%
        dplyr::summarise(tmm = sum(mm))

      out_df[i, 3] <- df3[[1]]

    }

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
#' @examples
#' \dontrun{
#' calc_mismatches(slist = data_list)
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
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
#' @examples
#' \dontrun{
#' dist <- dissimilarity(filtered_alleles = "numerical_alleles_filtered_a0.8_m0.2.csv",
#'     errors = "sample_error_results.csv")
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @import dplyr
#' @import tidyr
#' @import here
#' @importFrom readr read_csv
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
#' \code{misassign} compares different values of h and effect
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and a maximum number of allowable mismatches
#'     to generate a series of plots to visualise the two distributions of inclusion
#'     status to a group (in or out) based on number of loci at different mismatch
#'     thresholds.
#'
#'     Misassignment is defined as occurring when the 99.5% quantile of group
#'     inclusion (red dash line in plots) exceeds the 0.05% quantile of non-inclusion
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
#' @examples
#' \dontrun{
#'misassign(dist = dissimilarity_list, maxh = 5)
#' }
#'
#' @author Rujiporn Sun, \email{rujiporn.sun@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
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
#' @examples
#' \dontrun{
#'group_data <- group_membership(dist = dissimilarity_list, h = 5,
#'filtered_alleles = "numerical_alleles_filtered_a0.8_m0.2.csv")
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom readr read_csv write_csv
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
                                paste0("hclust_numerical_mismatch_",
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
#' @examples
#' \dontrun{
#'majorities(dist = dissimilarity_list, h = 5,
#'filtered_alleles = "numerical_alleles_filtered_a0.8_m0.2.csv",
#'errors = "sample_error_results.csv")
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://dbca-wa.github.io/DBCAscatR/index.html}
#' {the DBCAscatR website}
#'
#' @importFrom readr read_csv write_csv
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
                                paste0("hclust_numerical_mismatch_", h,
                                       "_group_majorities_and_ties.csv")))
  })

}
