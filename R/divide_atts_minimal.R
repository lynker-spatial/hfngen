#' Compute minimal divide attributes from rasters and GW parquet
#'
#' @description
#' Derives a minimal set of **divide-level attributes** by running zonal
#' statistics on input rasters and area-weighting a groundwater (Zmax) field
#' from a parquet dataset, then joins results by `divide_id`.
#'
#' @details
#' It computes:
#' - **Mode ISLTYP/IVGTYP** per divide using `zonal::execute_zonal()` on the
#'   provided rasters (via `terra::rast()`).
#' - **Mean impervious fraction** (`IMP`) per divide (zonal mean).
#' - **Area-weighted mean Zmax** per divide by joining GW parquet data
#'   (`ComID`, `Zmax`) to the `network` layer on
#'   `reference_id` = `ComID`, weighting by `incremental_areasqkm`.
#'
#' The results are returned as a tibble.
#'
#' @param gpkg Character. Path to a GeoPackage containing `divides` and `network` layers.
#' @param ISLTYP Character. Path to a raster of soil type (ISLTYP).
#'   Defaults to a CONUS path.
#' @param IVGTYP Character. Path to a raster of vegetation type (IVGTYP).
#'   Defaults to a CONUS path.
#' @param IMP Character. Path to a raster for impervious fraction used for
#'   zonal mean.
#' @param GW DuckDB connection via `hfutils::tbl_http()` with groundwater attributes,
#'   expected to include columns `ComID` and `Zmax`.
#'
#' @return A tibble with one row per `divide_id` containing:
#' - mode-based ISLTYP/IVGTYP attributes (column names as returned by `zonal`),
#' - `IMP` zonal mean,
#' - area-weighted mean `Zmax` (missing values imputed to 10).
#'
#' @section Requirements:
#' - The `divides` layer must contain `divide_id` and valid geometries.
#' - The `network` layer must contain:
#'   `divide_id`, `incremental_area_id`, `areasqkm`,
#'   `incremental_areasqkm`, and `reference_id`.
#' - Rasters should be in a compatible projection or reprojectable by `terra`.
#'
#' @section Imputation:
#' Any `Zmax` remaining `NA` after joins/weighting is imputed to **10**.
#'
#' @examples
#' \dontrun{
#' out <- minimal_div_attts(
#'   gpkg   = "path/to/hydrofabric.gpkg",
#'   ISLTYP = "path/to/ISLTYP.tif",
#'   IVGTYP = "path/to/IVGTYP.tif",
#'   IMP    = "path/to/IMP.tif",
#'   GW     = "path/to/nwm_gw_basins.parquet"
#' )
#' dplyr::glimpse(out)
#' }
#' @importFrom hfutils as_ogr st_as_sf
#' @importFrom terra rast
#' @importFrom dplyr select distinct collect left_join group_by summarize across mutate
#' @importFrom tibble as_tibble
#' @importFrom powerjoin power_full_join
#' @importFrom arrow open_dataset
#' @importFrom stats weighted.mean
#' @importFrom zonal execute_zonal
#' @export


minimal_div_atts = function(
    gpkg,
    ISLTYP = 'gridded/nwm/CONUS/ISLTYP.tif',
    IVGTYP = 'gridded/nwm/CONUS/IVGTYP.tif',
    IMP    = 'gridded/nwm/CONUS/derived/CONUS_imp_240.tif',
    GW     = 'tabular/nwm_gw_basins.parquet'){

  geom <- as_ogr(gpkg, 'divides') |>
    st_as_sf()

  nwm <-  suppressWarnings({
    zonal::execute_zonal(rast(c(ISLTYP, IVGTYP)),
                         geom,
                         fun = 'mode',
                         ID = "divide_id",
                         join = FALSE) |>
      setNames(c("divide_id", "mode.ISLTYP", "mode.IVGTYP"))
  })

  imp <-   suppressWarnings({
    zonal::execute_zonal(rast(IMP),
                         geom,
                         fun = 'mean',
                         ID = "divide_id",
                         join = FALSE) |>
      setNames(c("divide_id", "mean.impervious"))
  })

  zmax_data <-  GW |>
    select(ComID, Zmax) |>
    collect()

  zmax <-  as_ogr(gpkg, 'network') |>
    select(divide_id, incremental_area_id, areasqkm, incremental_areasqkm, reference_id) |>
    distinct() |>
    collect() |>
    left_join(zmax_data, by = c('reference_id' = 'ComID')) |>
    group_by(divide_id) |>
    summarize(across(c('Zmax'), ~ round(
      weighted.mean(.x, w = incremental_areasqkm, na.rm = TRUE), 5
    ))) |>
    mutate(mean.Zmax = ifelse(is.na(Zmax), 10, Zmax),
           Zmax = NULL)

  powerjoin::power_full_join(list(nwm, imp, zmax), by = 'divide_id') |>
    tibble::as_tibble()
}





