# Main plotting functions

#' Amplification plots
#'
#' \code{amp_splots} plots sample amplification rate against sample error rates.
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
#'@examples
#'\dontrun{
#'amp_splots()
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
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
      ggsave(fname, p3)
      print(p3)
    } else {
      stop(paste0("Not multivariate. Data probably has no error values and was ",
                  "most likely generated from raw data with no replicates"))
    }
  })

}

#' Data "missingness" plot
#'
#' \code{miss_hist} creates a histogram of the proportion missing data (NA) per locus.
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
#'@examples
#'\dontrun{
#'miss_hist()
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
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
      labs(title = "Locus missingness histogram",
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
#'@examples
#'\dontrun{
#'amp_hist()
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
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
    ggsave(here::here("results", "threshold", "amplification_histogram.jpg"), p1)
    print(p1)
  })
}

#' Dendrogram plot
#'
#' \code{dendro_plot} creates a dendrogram of samples clustered by genotype similarity.
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
#'@examples
#'\dontrun{
#'dendro_plot(dist = dissimilarity_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
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
#'@examples
#'\dontrun{
#'elbow_plot(dist = dissimilarity_list, maxh = 5)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
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
#' \code{freq_hist} creates a histogram of the frequency of genotype mismatches.
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
#'@examples
#'\dontrun{
#'freq_hist(dist = dissimilarity_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @importFrom tibble tibble
#' @import ggplot2
#' @import here
#'
#' @export
freq_hist <- function(dist){
  suppressWarnings({
    mmobj =  tibble::tibble(id = "s",
                            missmatch = dist[['combo']]$mm)

    p1 <- ggplot(aes(x= missmatch), data = mmobj) +
      geom_histogram(colour = "black", fill = "grey") +
      stat_bin(geom="text", colour="black", size=3.5,
               aes(label=..count..), position=position_stack(vjust=1.2)) +
      labs(y = "Frequency",
           x = "Missmatch",
           title = "Histogram of missmatch frequency") +
      theme_bw()
    ggsave(here::here("results", "cluster",
                      "hclust_group_mismatch_histogram.png"), p1)
    print(p1)

  })
}

#' Mismatch heatplot
#'
#' \code{heat_plot} creates a heat plot based on mismatches.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and creates a heat map. The heat map can
#'     help visualise groupings of samples based on mismatch values. It generates
#'     a "plotly" visualisation in the viewer when using RStudio whereby values
#'     can be interrogated by hovering the mouse cursor. The heat map can be saved
#'     to file either by using the interactive menu in the "plotly" html or by using
#'     the export button in the RStudio viewer. The html version is saved to
#'     `results/cluster/`.
#'
#' @inheritParams dendro_plot
#'
#' @return It will save and print to screen a plotly interactive html heat map.
#'
#'@examples
#'\dontrun{
#'heat_plot(dist = dissimilarity_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @importFrom heatmaply heatmaply
#' @import here
#'
#' @export
heat_plot <- function(dist){
  suppressWarnings({
    suppressMessages({
      distobj <- dist[['dist']]
      m <- as.matrix(distobj)
      heatmaply::heatmaply(m, k_row = 2, k_col = 2,
                           file = here::here("results", "cluster", "hclust_heatmap_plot.html"),
                           width = 1200, height = 750, plot_method = "ggplot")
    })
  })

}

#' Sample correlation heatplot
#'
#' \code{heat_cor_plot} creates a heat plot based on correlations.
#'
#' @details When run it uses the list object created by running
#'     \code{\link{dissimilarity}} and generates a correlation heat map. The
#'     correlation heat map can help visualise groupings of samples based on
#'     mismatch values. It generates a "plotly" visualisation in the viewer when
#'     using RStudio whereby values can be interrogated by hovering the mouse
#'     cursor. The heat map can be saved to png file either by using the interactive
#'     menu in the "plotly" html or by using the export button in the RStudio
#'     viewer. The html version is saved to `results/cluster/`.
#'
#' @inheritParams dendro_plot
#'
#' @return It will save and print to screen a plotly interactive html correlation
#'     heat map.
#'
#'@examples
#'\dontrun{
#'heat_cor_plot(dist = dissimilarity_list)
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @importFrom heatmaply heatmaply_cor
#' @import here
#'
#' @export
heat_cor_plot <- function(dist){
  suppressWarnings({
    suppressMessages({
      distobj <- dist[['dist']]
      m <- as.matrix(distobj)
      heatmaply::heatmaply_cor(m, k_row = 2, k_col = 2,
                               file = here::here("results", "cluster",
                                                 "hclust_heatmap_cor_plot.html"),
                               width = 1200, height = 750, plot_method = "ggplot")
    })

  })
}

#' Generate a leaflet plot showing site and individuals locations.
#'
#' \code{leaflet_map} generates a self contained leaflet html map showing
#'     locations of capture sites and individuals.
#'
#' @details When run it it creates an html leaflet map showing the locations of
#'     the field sites and where the individuals were sampled. The map is
#'     zoomable, layers can be turned on and off and there is a widget for
#'     calculating approximate distances and areas. Please note that the accuracy
#'     of the measurements is calculated using geodetic coordinates and is
#'     unlikely to be as accurate as measurements calculated using projected
#'     coordinates. The map will be saved to `results/finalised/`.
#'
#' @inheritParams summary_tables
#'
#' @return It will write to html a leaflet map showing site and individuals locations.
#'
#'@examples
#'\dontrun{
#'leaflet_map(groups_csv = "hclust_numerical_mismatch_h4_withGroups.csv",
#'metadata = "lookup.csv", prefix = "ID_", sample = "sample", site_ID = "roost_name",
#'field_date = "collection_date", lat = "dec_lat", long = "dec_long")
#'}
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' For more details see  \url{https://bartesto.github.io/ScatMatch/index.html}
#' {the ScatMatch website}
#'
#' @import readr
#' @importFrom stringr str_split
#' @importFrom sf st_as_sf
#' @importFrom htmlwidgets saveWidget
#' @import dplyr
#' @import lubridate
#' @import leaflet
#' @import here
#'
#' @export
leaflet_map <- function(groups_csv, metadata, prefix, sample, site_ID, field_date, lat, long){
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

    ## map
    # get coords for each site
    site_dat <- d2 %>%
      dplyr::select(site, lat, long) %>%
      dplyr::distinct() %>%
      sf::st_as_sf(coords = c("long", "lat"), crs = 4326)
    # get coords for each individual
    ind_dat <- d2 %>%
      dplyr::select(group, lat, long) %>%
      dplyr::rename(individual = group) %>%
      dplyr::distinct() %>%
      sf::st_as_sf(coords = c("long", "lat"), crs = 4326)

    map <- leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
      # setView(lat = -23, lng = 119, zoom = 11) %>%
      addProviderTiles(providers$Esri.WorldImagery) %>%
      addMiniMap(zoomLevelFixed = 2) %>%
      addMarkers(data = site_dat, popup = ~as.character(site),
                 label = ~as.character(site),
                 clusterOptions = markerClusterOptions(),
                 clusterId = "site",
                 group = "sites") %>%
      addMarkers(data = ind_dat, popup = ~as.character(individual),
                 label = ~as.character(individual),
                 clusterOptions = markerClusterOptions(),
                 clusterId = "individual",
                 group = "individuals") %>%
      addLayersControl(overlayGroups = c("sites", "individuals")) %>%
      addMeasure(primaryLengthUnit = "meters",
                 primaryAreaUnit = "hectares")

    htmlwidgets::saveWidget(map, file = here::here("results",
                                      "finalised",
                                      paste0("site_individual_locations_",
                                             ngrps, ".html")))

  })
}
