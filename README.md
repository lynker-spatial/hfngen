hfngen
================

`hfngen` provides utilities for turning Hydrofabric GeoPackages into
NextGen-ready datasets. It focuses on three workflows:

1.  `build_nextgen()` – convert Hydrofabric flowpath/flowline packages
    into a prefixed NextGen network with optional flowline detail and
    hydrolocation enrichment.
2.  `minimal_div_atts()` – derive a minimal set of divide-level
    attributes by combining zonal raster statistics (soil, vegetation,
    imperviousness) with groundwater parquet data.
3.  `calculate_flow_atts()` – roll up feature-level hydraulic attributes
    (served through Arrow/DuckDB) to the flow identifier of your choice
    with length-weighted means.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("lynker-spatial/hfngen")
```

## Data expectations

- **Hydrofabric inputs** must include `flowpaths`, `divides`, and
  (optionally) flowline-resolution layers (`flowlines`,
  `incremental_areas`, `hydrolocations`).
- **Raster inputs** should share a compatible CRS with the divides
  layer. The package will reproject geometries on the fly when needed,
  but aligning rasters ahead of time avoids repeated transforms.
- **Hydraulic attributes** need a `feature_id` column that matches
  `network$reference_id`, the attributes you wish to aggregate, and a
  length-based weighting column.

## Support

Please open an issue or pull request if you run into problems or would
like to contribute additional workflows.
