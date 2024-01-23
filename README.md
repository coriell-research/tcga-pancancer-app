## TCGA PanCancer Data Explorer

A simple Shiny app to explore data from the TCGA PanCancer Atlas.

## Data Source

Data used in this app was downloaded from the [TCGA PanCancer Atlas Hub](https://xenabrowser.net/datapages/?host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) via UCSC Xena browser. The `scripts/01_download.sh` file provides paths. Data was then parsed into a format amenable to exploration with a relational database (duckdb) using the `scripts/02_create-bd.qmd` script. This data is then exposed to the user with a simple Shiny app interface.

