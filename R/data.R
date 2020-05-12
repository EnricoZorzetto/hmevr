#' Series of daily precipitation recorded in New York City Central Park
#' part of the NOAA USHCN dataset
#'
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{PRCP}{Daily precipitation, in mm/day}
#'   \item{DATE}{Date in integer format YYYYMMDD}
#'   \item{YEAR}{Year in integer format YYYY}
#'   \item{MONTH}{Month in format MM}
#'   \item{DAY}{Day in format DD}
#'   \item{STATION}{Station indentification number}
#'   \item{QFLAG}{Quality flag}
#'   \item{QFLAG}{Missing data flag}
#'   \item{SFLAG}{Source flag}
#'   ...
#' }
#' @source \url{https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/us-historical-climatology-network-ushcn}
"nycp"
