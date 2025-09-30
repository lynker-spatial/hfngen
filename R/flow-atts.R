#' Calculate flow-level hydraulic attributes by length-weighted aggregation
#'
#' Aggregates channel/hydraulic attributes from a feature-level attributes dataset
#' (served via **Arrow**) up to a hydrologic flow identifier (e.g., `flowline_id`)
#' using a length-based weighting column. The result returns common hydraulic
#' terms (bottom/top width, side slope, composite roughness, etc.) along with any
#' gage references attached to the flow.
#'
#' @param network A data frame/tibble containing at least:
#'   - `reference_id` (IDs that link to the attribute dataset),
#'   - `hl_class` and `hl_reference` (to extract gage references),
#'   - the grouping column given by `hf_id` (default: `"flowline_id"`),
#'   - and (optionally) other columns used downstream.
#'   This is typically the `network` layer produced by your NextGen build.
#'
#' @param hydraulics_path DuckDB connection via `hfutils::tbl_http()`
#'   - `feature_id` (will be mapped to `reference_id`),
#'   - the columns listed in `attribute_cols`,
#'   - the weighting column named in `weight_col`.
#'
#' @param hf_id Character name of the network flow identifier column to group by.
#'   Defaults to `"flowline_id"`. Must be present in `network`.
#'
#' @param attribute_cols Character vector of attribute column names to aggregate
#'   (present in the Arrow dataset). Defaults to
#'   `c('owp_y_bf','owp_tw_bf','bf_area','owp_roughness_bathy','slope')`.
#'   These are averaged with a length-weighted mean.
#'
#' @param weight_col Character name of the length-based weighting column present
#'   in the Arrow dataset. Defaults to `"lengthkm"`.
#'
#' @return A tibble with one row per `hf_id` containing:
#'   - `{hf_id}`: the flow identifier
#'   - `BtmWdth`: estimated bottom width (m)
#'   - `TopWdth`: estimated top width (m)
#'   - `ChSlp`:  side slope (horizontal:vertical), floored at `1e-5`
#'   - `nCC`:    composite roughness (simple 2× `owp_roughness_bathy`)
#'   - `So`:     channel slope (from `slope`)
#'   - `TopWdthCC`: composite-channel top width (3× `TopWdth`)
#'   - `Y`:      bankfull depth (from `owp_y_bf`)
#'   - `hl_reference`: comma-separated gage references associated to the flow
#'
#' All weighted means ignore `NA`s (`na.rm = TRUE`). Be sure the weighting column
#' (`weight_col`) reflects segment length in consistent units with your study.
#'
#' @examples
#' \dontrun{
#' # Arrow dataset (parquet directory) with feature_id, lengthkm, and OWP attributes
#' attr_path <- "s3://bucket/owp/channel_atts_parquet/"
#'
#' # Network table produced earlier, with columns:
#' # reference_id, hl_class, hl_reference, flowline_id, ...
#' net <- readr::read_csv("network.csv")
#'
#' atts <- calculate_flow_atts(
#'   network        = net,
#'   attribute_path = attr_path,
#'   hf_id          = "flowline_id",
#'   attribute_cols = c("owp_y_bf","owp_tw_bf","bf_area","owp_roughness_bathy","slope"),
#'   weight_col     = "lengthkm"
#' )
#' }
#'
#' @importFrom dplyr select filter collect mutate distinct left_join group_by summarise across
#' @importFrom arrow open_dataset
#' @importFrom stats weighted.mean
#' @export

calculate_flow_atts <- function(network,
                                hydraulics,
                                hf_id = "flowline_id",
                                attribute_cols = c('owp_y_bf', 'owp_tw_bf','bf_area',
                                                   'owp_roughness_bathy','slope'),
                                weight_col = 'lengthkm') {

  map <- hydraulics |>
    dplyr::select(reference_id = feature_id,
                  all_of(attribute_cols),
                  small_w = all_of(weight_col)) |>
    dplyr::filter(reference_id %in%  !!unique(network$reference_id)) |>
    dplyr::collect()

 idx <- dplyr::select(network, reference_id, hl_class, hl_reference, all_of(hf_id)) |>
    dplyr::mutate(hl_reference = ifelse(hl_class == 'gage', hl_reference, NA)) |>
    dplyr::select(-hl_class) |>
    dplyr::distinct() |>
    dplyr::left_join(map, by = 'reference_id')  |>
    dplyr::group_by(!!dplyr::sym(hf_id)) |>
    dplyr::mutate(hl_reference = paste(na.omit(hl_reference), collapse = ',')) |>
    dplyr::summarise(dplyr::across(
      .cols = all_of(attribute_cols),
      .fns = ~ weighted.mean(.x, w = small_w, na.rm = TRUE),
      .names = "{.col}"
    ),
    hl_reference = hl_reference[1],
    .groups = 'drop') |>
    dplyr::mutate(
      Bw = (2 * bf_area / owp_y_bf) - owp_tw_bf,
      Bw = ifelse(Bw <= 0, .1, Bw)
    ) |>
    dplyr::select(
      all_of(hf_id),
      n = owp_roughness_bathy,
      Tw = owp_tw_bf,
      Bw,
      Y = owp_y_bf,
      So = slope,
      hl_reference
    ) |>
    dplyr::mutate(
      BtmWdth = pmin(Bw, Tw),
      TopWdth = pmax(Bw, Tw),
      nCC = 2 * n,
      TopWdthCC = 3 * TopWdth,
      ChSlp = ((TopWdth - BtmWdth) / 2) / Y,
      ChSlp = ifelse(ChSlp <= 1e-5, 1e-5, ChSlp),
      BtmWdth = ifelse(ChSlp == 1e-5, TopWdth - (2 * ChSlp * Y), BtmWdth),
      MusX = 0.2,
      MusK = 3600.0
    ) |>
    dplyr::select(hf_id,
           BtmWdth,
           TopWdth,
           ChSlp,
           n,
           nCC,
           So,
           TopWdthCC,
           hl_reference,
           MusX, MusK)

 dplyr::select(network, dplyr::all_of(hf_id), lengthkm = !!gsub("id", "lengthkm", hf_id)) |>
   dplyr::distinct() |>
   dplyr::mutate(lengthm = lengthkm * 1000, lengthkm = NULL) |>
   dplyr::right_join(idx, by = hf_id)


}

