#' Build a NextGen-style network from a Hydrofabric GeoPackage
#'
#' Translates a Hydrofabric flowpath/divide GeoPackage into a NextGen-style
#' topology with prefixed identifiers, optional flowline-level expansion, and an
#' integrated `network` table suitable for downstream modeling and ETL.
#'
#' @param flowpath_gpkg Character path to the **flowpath/divide** Hydrofabric
#'   GeoPackage. Must include at least layers `flowpaths` and `divides`.
#' @param flowline_gpkg  Optional character path to a **flowline-resolution**
#'   Hydrofabric GeoPackage. If provided, `flowlines` and `incremental_areas`
#'   tables are created and the `network` table is enriched with hydrolocation
#'   attributes where available.
#' @param collapse_method Character, one of `"upstream"` (default) or
#'   `"downstream"`. Controls which side receives the nexus identifier when
#'   collapsing node-to-edge relationships.
#' @param nexus_prefix Character prefix for created nexus identifiers.
#' @param terminal_nexus_prefix Character prefix for created **terminal** nexus
#'   identifiers. (Reserved; currently all nexuses use `nexus_prefix`.)
#' @param coastal_nexus_prefix Character prefix for created **coastal** nexus
#'   identifiers. (Reserved.)
#' @param internal_nexus_prefix Character prefix for created **internal** nexus
#'   identifiers. (Reserved.)
#' @param divide_prefix Character prefix for divide (catchment) identifiers.
#' @param flowpath_prefix Character prefix for flowpath identifiers.
#' @param term_add Numeric scalar added to terminal `flowpath_toid` to create a
#'   unique terminal-nexus keyspace (defaults to 1e9).
#' @param outfile Optional character file path (\*.gpkg). When provided, writes
#'   all resulting layers to the specified GeoPackage; otherwise returns an R
#'   list with the layers.
#'
#' @return Invisibly returns either the written `outfile` path (when `outfile`
#'   is provided) or a named list with layers (`nexus`, `flowpaths`, `divides`,
#'   `flowlines`*, `incremental_areas`*, `network`), where `*` indicates layers
#'   present only when `flowline_gpkg` is supplied.
#' @examples
#' \dontrun{
#' out <- build_nextgen(
#'   flowpath_gpkg = "hf_ref.gpkg",
#'   flowline_gpkg  = "hf_flowline.gpkg",
#'   outfile        = "nextgen.gpkg"
#' )
#' }
#'
#' @importFrom sf st_drop_geometry st_as_sf st_coordinates st_intersects
#' @importFrom dplyr select filter mutate arrange distinct group_by ungroup left_join right_join bind_cols bind_rows if_else full_join
#' @importFrom tidyr separate_longer_delim
#' @importFrom cli cli_alert_success cli_alert
#' @importFrom hfutils read_hydrofabric write_hydrofabric layer_exists get_node node_geometry
#' @export


build_nextgen = function(flowpath_gpkg,
                         flowline_gpkg = NULL,
                         collapse_method       = "upstream",
                         nexus_prefix          = "nex-",
                         terminal_nexus_prefix = "tnx-",
                         coastal_nexus_prefix  = "cnx-",
                         internal_nexus_prefix = "inx-",
                         divide_prefix   = "cat-",
                         flowpath_prefix = "fp-",
                         term_add        = 1e9,
                         outfile = NULL){

  cli::cli_alert("\n--- Applying NextGen topology ---\n")

  fp_network_list  <- read_hydrofabric(flowpath_gpkg)

  ngen_network_list <- list()

  fp_net       <- select(fp_network_list$flowpaths,
                         flowpath_id,
                         flowpath_toid,
                         flowline_id,
                         areasqkm,
                         hydroseq)


  fp_term_net  <- filter(fp_net, flowpath_toid == 0 | flowpath_toid >= term_add)

  terms <- filter(fp_net, flowpath_id %in% fp_term_net$flowpath_id)

  term_net <- bind_cols(
    terms,
    st_coordinates(get_node(terms, "end"))
  ) |>
    group_by(X, Y) |>
    mutate(flowpath_toid = dplyr::cur_group_id() + term_add) |>
    ungroup() |>
    mutate(X = NULL, Y = NULL)

  nodes  <- node_geometry(fp_net)

  nodes$tofrom <- lengths(st_intersects(nodes, node_geometry(fp_net, position = "start")))  > 0

  nodes$hw <- !nodes$flowpath_id %in% nodes$flowpath_toid

  imap  <- st_intersects(fp_net, nodes)

  topo  <- data.frame(
    flowpath_id    = nodes$flowpath_id[unlist(imap)],
    flowpath_toid  = rep(fp_net$flowpath_id, times = lengths(imap)),
    hs             = nodes$hydroseq[unlist(imap)],
    tofrom         = nodes$tofrom[unlist(imap)],
    hw         = nodes$hw[unlist(imap)]) |>
    filter(flowpath_id != flowpath_toid) |>
    group_by(flowpath_toid) |>
    mutate(flag = any(hw == TRUE)) |>
    ungroup()

  hw_nodes = filter(topo, flag == TRUE) |>
    dplyr::arrange(flowpath_toid, -hs) |>
    distinct(flowpath_toid, .keep_all = TRUE)

  topo = filter(topo, flag != TRUE) |>
    dplyr::arrange(flowpath_toid,
                   desc(tofrom),
                   desc(hs)) |>
    distinct(flowpath_toid, .keep_all = TRUE) |>
    bind_rows(select(st_drop_geometry(term_net), flowpath_id, flowpath_toid)) |>
    bind_rows(hw_nodes)

  if(collapse_method == "upstream"){
    ngen_network_list$nexus <-  dplyr::semi_join(nodes, topo, by = "flowpath_id") |>
      mutate(nexus_id   = if_else(is.na(flowpath_toid), NA, paste0(nexus_prefix,    flowpath_toid)),
             nexus_toid = if_else(is.na(flowpath_toid), NA, paste0(flowpath_prefix, flowpath_toid))) |>
      select(nexus_id, nexus_toid)
  } else {
    ngen_network_list$nexus <-  dplyr::semi_join(nodes, topo, by = "flowpath_id") |>
      mutate(nexus_id   = if_else(is.na(flowpath_id),   NA, paste0(nexus_prefix,    flowpath_id)),
             nexus_toid = if_else(is.na(flowpath_toid), NA, paste0(flowpath_prefix, flowpath_toid))) |>
      select(nexus_id, nexus_toid)
  }

  ngen_network_list$flowpaths <-
    fp_network_list$flowpaths |>
    mutate(flowpath_id   = if_else(is.na(flowpath_id),   NA, paste0(flowpath_prefix, flowpath_id)),
           flowpath_toid = if_else(is.na(flowpath_toid), NA, paste0(nexus_prefix, flowpath_toid)),
           divide_id = if_else(is.na(divide_id), NA, paste0(divide_prefix, divide_id)))

  ngen_network_list$divides <-
    fp_network_list$divides |>
    mutate(divide_id   = if_else(is.na(divide_id), NA, paste0(divide_prefix, divide_id)),
           flowpath_toid = NULL) |>
    left_join(st_drop_geometry(  ngen_network_list$flowpaths) |> select(divide_id, flowpath_id), by = "divide_id")

  stopifnot(
    hydrofab:::network_is_dag(ngen_network_list$flowpaths, "flowpath_id", "flowpath_toid")
  )

  if(!is.null(flowline_gpkg)){
    cli::cli_alert_success("\n--- Flowline topology provided ---\n")

    fl_network_list <- read_hydrofabric(flowline_gpkg)

    stopifnot(
      hydrofab:::network_is_dag(fl_network_list$flowpaths, "flowpath_id", "flowpath_toid")
    )


    ngen_network_list$flowlines = select(st_drop_geometry(ngen_network_list$flowpaths), flowpath_id, flowpath_toid, flowline_id, areasqkm) |>
      tidyr::separate_longer_delim(cols = 'flowline_id', delim = ",") |>
      mutate(flowline_id = as.integer(flowline_id)) |>
      left_join(select(st_drop_geometry(fl_network_list$flowpaths),
                       flowline_id = flowpath_id,
                       flowline_toid = flowpath_toid,
                       flowline_hydroseq = hydroseq,
                       area_incr = areasqkm),
                by = "flowline_id",
                relationship = "many-to-many") |>
      group_by(flowpath_id) |>
      mutate(flowline_incr = rank(flowline_hydroseq),
             area_ratio = round(area_incr / areasqkm, 3)) |>
      ungroup() |>
      select(flowline_id, flowline_toid,
             flowline_incr, flowline_hydroseq,
             areasqkm = area_incr, area_ratio,
             flowpath_id, flowpath_toid) |>
    left_join(select(fl_network_list$flowpaths,
                     flowline_id = flowpath_id,
                     incremental_area_id = divide_id,
                     reference_id = member_comid,
                     vpuid),
              by  = "flowline_id",
              relationship = "many-to-many") |>
      st_as_sf()

    ngen_network_list$incremental_areas <-
      select(st_drop_geometry(ngen_network_list$divides), flowpath_id, divide_id) |>
      dplyr::right_join(select(st_drop_geometry(ngen_network_list$flowlines), flowpath_id, incremental_area_id = flowline_id, area_ratio), by = "flowpath_id") |>
      left_join(select(fl_network_list$divides, incremental_area_id = divide_id, vpuid, areasqkm), by  = "incremental_area_id") |>
      st_as_sf()

  }

  if(hfutils::layer_exists(flowpath_gpkg, "pois")){
    ngen_network_list$pois = sf::read_sf(flowpath_gpkg, "pois") |>
      select(poi_id, flowpath_id) |>
      mutate(
        flowpath_id = if_else(is.na(flowpath_id), NA, paste0(flowpath_prefix, flowpath_id))
      )
  }

    net <- sf::read_sf(flowpath_gpkg, "network") |>
      select(flowpath_id, reference_id = hf_id) |>
      mutate(
        flowpath_id = if_else(is.na(flowpath_id), NA, paste0(flowpath_prefix, flowpath_id)),
      )

    net = ngen_network_list$flowpaths |>
      select(flowpath_id, flowpath_toid, hydroseq, divide_id, flowline_id, poi_id, flowpath_lengthkm = lengthkm, areasqkm, tot_drainage_areasqkm) |>
      st_drop_geometry() |>
      tidyr::separate_longer_delim(cols = 'flowline_id', delim = ",") |>
      mutate(flowline_id = as.integer(flowline_id)) |>
      dplyr::right_join(net, by = "flowpath_id", relationship = "many-to-many")

    if(!is.null(flowline_gpkg)){
      net = ngen_network_list$flowlines %>%
        mutate(flowline_lengthkm = hfutils::add_lengthkm(.)) |>
        select(flowline_id, flowline_toid, flowline_hydroseq, flowline_lengthkm, incremental_area_id, incremental_areasqkm = areasqkm) |>
        st_drop_geometry() |>
        dplyr::right_join(net, by = "flowline_id", relationship = "many-to-many")

    if(hfutils::layer_exists(flowline_gpkg, "hydrolocations")){

      hl <- sf::read_sf(flowline_gpkg, 'hydrolocations') |>
        select(any_of(c('poi_id', 'hl_reference', 'hl_class', 'hl_source', 'hl_uri'))) |>
        distinct() |>
        mutate(poi_id = as.character(poi_id))

      net <- dplyr::right_join(st_drop_geometry(hl), net, by = "poi_id", relationship = "many-to-many")

      ngen_network_list$hydrolocations <- hl
    }

    ngen_network_list$network = select(net,
           starts_with("flowpath"),
           'hydroseq',
           "divide_id",
           starts_with("flowline"),
           "incremental_area_id",
           'incremental_areasqkm',
           "poi_id",
           starts_with("hl"),
           'reference_id',
           dplyr::everything()
           ) |>
      distinct()
    }

    if(!is.null(outfile)){
      write_hydrofabric(ngen_network_list, outfile = outfile, enforce_dm = FALSE)
      if("hl" %in% names(ngen_network_list)){
        cli::cli_alert_success("Adding hydrolocations written to {.file {outfile}}")
        write_sf(ngen_network_list$hydrolocations, outfile, layer = "hydrolocations")
      }
    } else {
      return(invisible(ngen_network_list))
    }

}
