#-------------------------------------------------------------------------------
#
#     Accumulation of local biodiversity information across the tree of life 
#       (using ERCCIS data as an example) - Revision 1
#         By David Baker (11/2021)
#
#-------------------------------------------------------------------------------

#' This script analyses the accumulation of local biodiversity records at fine 
#' spatial resolutions (finest scale 1km2) for 15 species groups over the period 
#' 1960 to 2019 and evaluates the spatial drivers (socio-environmental) of 
#' sampling completeness.

#-------------------------------------------------------------------------------
#
# SECTION 1. Setup
#
#-------------------------------------------------------------------------------

#-- Clear environment
rm(list = ls())

#-- Set to UTC
Sys.setenv(TZ = "UTC")

#-- Increase memory size
memory.size(1000000)

#-- Set seed
set.seed(17)

#-- Path to source functions
function_path <- "Scripts/Functions_main_R1_fs.R"
source(function_path)

#-- Shared data path
shareDir <-
  "C:/Users/db625/OneDrive - University of Exeter/Fellowship_data/"

#-- Load packages
packageList <- c(
  "snowfall",
  # Spatial data package
  "rgdal",
  "raster",
  "rgeos",
  "sf",
  "spdep",
  # Parsing date and time 
  "lubridate",
  # Data handling
  "tidyr",
  "tidylog",
  "data.table",
  "RSQLite",
  "dplyr",
  "dtplyr",
  "gtools",
  "ncf",
  # plotting
  "ggplot2",
  "grid",
  "gridExtra",
  "cowplot",
  "RColorBrewer",
  "viridis",
  "scico",
  "ggthemes",
  # Statistical modelling
  "spatstat",
  "brms",
  "rareNMtests",
  "Rarefy",
  "mgcv",
  "gratia",
  "DHARMa",
  "iNEXT",
  "iNEXT.4steps"
)
package.load(packageList)

#-- UK shapefile
cornwall <- readOGR(paste0(shareDir,
                           "EnglishCountiesCornwall_2011"),
                    "england_ct_2011")
cornwall_sf <- st_as_sf(cornwall)

#-- Species groups
sppNameTab <- tribble (
  
  ~SpeciesGroup, ~SpeciesGroupFile,
  #--|--|--
  "fungus",                      "Fungi",
  "lichen",                      "Lichens",
  "moss",                        "Mosses",
  "fern",                        "Ferns",
  "flowering plant",             "Flowering_Plants",
  "terrestrial mammal",           "Terrestrial_Mammals",
  "terrestrial mammal - bat (Chiroptera)", "Bats"  ,
  "amphibian",                   "Herptiles",
  "reptile",                     "Herptiles",
  "bird",                        "Birds", 
  "bony fish (Actinopterygii)",  "Bony_Fish",
  "insect - moth",                "Moths",
  "insect - butterfly",           "Butterflies",
  "insect - hymenopteran",        "Hymenoptera",
  "insect - true fly (Diptera)",  "Diptera",
  "insect - dragonfly (Odonata)", "Dragonflies"
  
)
sppNameTab <- as.data.frame(sppNameTab)

#-- Hexagonal grids
hexGrain1 <- makeHexGrid(cornwall_sf, 1)
hexGrain4 <- makeHexGrid(cornwall_sf, 4)
hexGrain16 <- makeHexGrid(cornwall_sf, 16)
hexGrain64 <- makeHexGrid(cornwall_sf, 64)

#-- Things to always keep in the global environment
keep_files <- c("keep_files",
                "shareDir",
                "cornwall_sf",
                "sppNameTab",
                "source.keep",
                "function_path",
                "hexGrain1",
                "hexGrain4",
                "hexGrain16",
                "hexGrain64")

#-------------------------------------------------------------------------------
#
# SECTION 2. Accumulation of records and the median age of records
#
#-------------------------------------------------------------------------------

#' We want to know how records accumulate over time so that we can see the 
#' temporal profile across space and time. Records in some areas might be older 
#' on average the in other, so are we comparing like-for-like. What drives these
#' differences spatially?

#---- Source functions for this section
source.keep(
  c(
    "erccis.db.tbl",
    "db.transform.proj",
    "species.group.select",
    "makeHexGrid",
    "recToHexSetUp",
    "recXgridXtime",
    "firstDeriv_nRec"
  ),
  function_path
)

#---- Load ERCCIS SQLite DB and tbl
erccis_data <- erccis.db.tbl(
  sqlite_file = paste0(shareDir, "erccis_oct2020.sqlite"),
  table_name = "erccis_oct2020",
  sppGrps = sppNameTab$SpeciesGroup,
  sppNameLookUp = sppNameTab
)

#---- Distribution of recording intensity over time and space
#-- Get number of records per species group, year, and grain
for (sppGrpN in unique(sppNameTab[, 2])) {
  print(sppGrpN)
  
  #~# Subset to just species group
  sppSelGrp <- NULL
  spp_grp <-
    as.character(sppNameTab[, 1][which(sppNameTab[, 2] == sppGrpN)])
  erccis_sppGrp <- species.group.select(erccis_data,
                                        spp_grp,
                                        cornwall_sf,
                                        "Cornwall",
                                        sppSelect = sppSelGrp)
  
  #~# Extract records per grain grn = 64
  for (grn in c(1, 4, 16, 64)) {
    hexGrain <- get(paste0(hexGrain, grn))
    NRecSpp <-
      recXgridXtime(
        erccis_sppGrp,
        sppGrp = sppGrpN,
        timeSlice = 1,
        region_sf = cornwall_sf,
        hexCellArea = grn,
        bioHex = NULL,
        hexGrid = hexGrain
      )
    save(
      NRecSpp,
      file = paste0(
        "Outputs/Number_Records_SppGrp/NRecSpp_",
        sppGrpN,
        "_grain_",
        grn,
        ".Rdata"
      )
    )
    
  }
  
}

#-- Extract to dataframe
listPth <- list.files("Outputs/Number_Records_SppGrp/")
sppGrp_N_Rec <- lapply(listPth, function(sppGrp) {
  sppDat <-
    get(load(paste0(
      "Outputs/Number_Records_SppGrp/", sppGrp
    )))
  sppDat$grain <-
    as.numeric(last(strsplit(sub(
      '.Rdata', '', sppGrp
    ),
    split = '_|\\.')[[1]]))
  sppDat
})
sppGrp_N_Rec <- do.call(bind_rows, sppGrp_N_Rec)
sppGrp_N_Rec$Yr <- gsub('.{5}$', '', sppGrp_N_Rec$TS)
sppGrp_N_Rec$sppGrp <- sub('_', ' ', sppGrp_N_Rec$sppGrp)

#-- Regional record number trend
#- Get yearly totals 
reg_n_trend <- sppGrp_N_Rec %>%
  filter(grain == 64) %>%
  select(!c("n")) %>%
  drop_na() %>%
  group_by(Yr,sppGrp) %>%
  summarise(nRec = sum(N), .groups = 'drop') %>%
  group_by(sppGrp) %>%
  mutate(cumRec = cumsum(nRec)) %>%
  arrange(sppGrp, Yr)

#- Get median year
reg_n_med <- reg_n_trend %>%
  group_by(sppGrp) %>%
  filter(abs(cumRec - max(cumRec) / 2) ==
           min(abs(cumRec - max(cumRec) / 2))) %>%
  arrange(sppGrp) %>%
  select(sppGrp, Yr) %>%
  rename(medRecYr = Yr)
modal(as.numeric(reg_n_med$medRecYr)) # modal year across species groups

#- Get peak year
peak_year <- reg_n_trend %>%
  group_by(sppGrp) %>%
  filter(nRec == max(nRec)) %>%
  arrange(sppGrp) %>%
  select(sppGrp, Yr, nRec) %>%
  rename(c(nRecMaxRecYr = nRec, MaxRecYr = Yr))

#- Get yearly totals for 1960
sppGrp_N_Rec_1960 <- reg_n_trend %>%
  filter(Yr == 1960) %>%
  mutate(cumRec = cumsum(nRec)) %>%
  arrange(sppGrp, Yr) %>%
  select(sppGrp, nRec) %>%
  rename(nRec1960 = nRec)

#- Get total 1960-2019
sppGrp_N_Rec_total <- reg_n_trend %>%
  filter(Yr == 2019) %>%
  arrange(sppGrp, Yr) %>%
  select(sppGrp, cumRec) %>%
  rename(totalRec = cumRec)

#- Get yearly totals
prop_cell_rec_2019 <- sppGrp_N_Rec %>%
  filter(grain %in% c(1, 4, 16, 64)) %>%
  arrange(sppGrp, id) %>%
  group_by(grain) %>%
  mutate(nCells = length(unique(id))) %>%
  group_by(sppGrp, grain, id, nCells) %>%
  summarise(nRec = sum(N), .groups = 'drop') %>%
  filter(nRec > 10) %>%
  group_by(sppGrp, grain, nCells) %>%
  summarise(cRec = n(), .groups = 'drop') %>%
  mutate(pc10 = round((cRec / nCells) * 100)) %>%
  select(sppGrp, grain, pc10) %>%
  pivot_wider(
    id_cols = sppGrp,
    names_from = grain,
    names_prefix = "grn",
    values_from = pc10
  )

#-- Create summary table
summaryTable <-
  Reduce(
    function(...)
      merge(..., by = 'sppGrp', all = TRUE),
    list(
      sppGrp_N_Rec_1960,
      sppGrp_N_Rec_total,
      reg_n_med,
      peak_year,
      prop_cell_rec_2019
    )
  )
summaryTable$sppGrp <- factor(summaryTable$sppGrp,
                              levels = sub("_", " ", unique(sppNameTab[, 2])))
summaryTable <- arrange(summaryTable, sppGrp)
mean(summaryTable$nRec1960, na.rm = TRUE)
mean(summaryTable$totalRec, na.rm = TRUE)
write.csv(summaryTable, 'Outputs/Figures_R1/Table_1_Summary.csv')

#-- First derivatives to identify the fastest acceleration
firstDerivRecNo <- lapply(unique(sppNameTab[, 2]), function(spp) {
  fd <- firstDeriv_nRec(reg_n_trend, sub('_', ' ', spp))
})
firstDerivRecNo <- do.call(bind_rows, firstDerivRecNo)
firstDerivRecNo %>% group_by(sppGrp) %>% summarise(yrM = min(Yr))

#-- Order factors
firstDerivRecNo$sppGrp <- 
  factor(firstDerivRecNo$sppGrp,
         levels = sub("_", " ", unique(sppNameTab[, 2])))
reg_n_trend$sppGrp <- 
  factor(reg_n_trend$sppGrp,
         levels = sub("_", " ", unique(sppNameTab[, 2])))
reg_n_med$sppGrp <- 
  factor(reg_n_med$sppGrp,
         levels = sub("_", " ", unique(sppNameTab[, 2])))

#----
#---- Plot for number of records per year for all species groups
#----
sppGrp_N_Rec_plot <- ggplot(reg_n_trend) +
  geom_segment(
    aes(
      x = as.numeric(medRecYr),
      y = -Inf,
      xend = as.numeric(medRecYr),
      yend = Inf
    ),
    colour = "#1F78B4",
    size = 1.5,
    linetype = 2,
    data = reg_n_med
  ) +
  geom_line(aes(
    x = as.numeric(Yr),
    y = log10(nRec),
    group = sppGrp
  ), linetype = 2) +
  geom_point(aes(x = as.numeric(Yr),  y = log10(nRec))) +
  geom_ribbon(aes(
    x = Yr,
    ymin = (lwr),
    ymax = (upr)
  ),
  alpha = 0.25, data = firstDerivRecNo) +
  geom_line(aes(x = Yr, y = (response)), data = firstDerivRecNo, size = 0.5) +
  geom_line(
    aes(
      x = Yr,
      y = (response),
      group = fdYrs
    ),
    colour = 'Red',
    size = 1,
    data = firstDerivRecNo[which(firstDerivRecNo$fdYrs > 0), ]
  ) +
  ylab(bquote("Number of records (log"['10'] * ")")) + xlab('Year') +
  theme_bw(base_size = 16) %+replace%
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(breaks = seq(1960, 2020, 20)) + ylim(0, 5) +
  facet_wrap( ~ sppGrp, ncol = 5)
ggsave(filename = "Outputs/Figures_R1/Figure_1_Record_density_time.png", 
       plot = sppGrp_N_Rec_plot, width = 12, height = 7)

#---- Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

#-------------------------------------------------------------------------------
#
# SECTION 3. Temporal accumulation of coverage: evaluate thresholds
#
#-------------------------------------------------------------------------------

#' We want to assess sampling completeness for these records and to do this 
#' we'll follow previous studies in using rarefaction curves to estimate 
#' sampling completeness. However, these estimates are likely to be sensitive
#' to the number of samples used to estimate completeness. We will, therefore,
#' use a simulation of the record collection methods to test the effect of 
#' sample size and assemblage size on sampling completeness estimates. The 
#' simulation models are derived from the fsSDM package.

#---- Load some environmental data for the region
#-- Load CEH land use data at 1 km2 with buffer 0 and 1
land_cover_c19 <-
  get(load(paste0(
    shareDir, "/cornwall_ceh2019_landcover_1km.Rdata"
  )))
#-- Load bioclimate variables
bioClim <-
  get(load(paste0(
    shareDir, "/cornwall_climvars_1km_means.Rdata"
  )))
#-- Merge into single enviro data frame
enviro_dat <- left_join(land_cover_c19, bioClim,
                        by = c('id', 'X', 'Y', 'grain'))
#-- Vector of environmental variables 
enviro_vars <- c("LC_1_0", "LC_3_0", "LC_4_0", "LC_7_0", "LC_9_0", "LC_10_0", 
                 "LC_11_0", "LC_14_0", "LC_15_0","LC_16_0", "LC_17_0", 
                 "LC_18_0", "LC_19_0", "LC_20_0", "LC_21_0", "temp_ann_m",
                 "temp_gs_m", "mtcm_m", "mtcq_m", "frost_m", "precip_gs_m",
                 "soilm_gs_m")

#---- Set up experiment
exp_runs <- expand.grid(
  species_n = c(10, 25, 50, 100, 200, 400),
  site_prop = c(0.05, 0.1, 0.2, 0.3),
  reliability = c('low', 'medium')
)
#-- Run exp
#- 4km - looks like >25 samples
exp_rarefy_4 <-
  test_coverage_est(exp_runs, 4, hex1)
write.csv(exp_rarefy_4, file = "Outputs/Test_SAC_4km.csv")
plot.sim.res("Outputs/Test_SAC_4km.csv", "4") 

#- 16km - looks like >50 samples
exp_rarefy_16 <-
  test_coverage_est(exp_runs, 16, hex1)
write.csv(exp_rarefy_16, file = "Outputs/Test_SAC_16km.csv")
plot.sim.res("Outputs/Test_SAC_16km.csv", "16")

#- 64km - looks like >50 samples
exp_rarefy_64 <-
  test_coverage_est(exp_runs, 64, hex1)
write.csv(exp_rarefy_64, file = "Outputs/Test_SAC_64km.csv")
plot.sim.res("Outputs/Test_SAC_64km.csv", "64")

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

#-------------------------------------------------------------------------------
#
# SECTION 4. Temporal accumulation of coverage
#
#-------------------------------------------------------------------------------

source.keep(
  c(
    "erccis.db.tbl",
    "db.transform.proj",
    "recToHexSetUp",
    "recXgridXtime",
    "siteXSpeciesXTime",
    "species.group.select",
    "run_cover",
    "ci_panel_plot",
    "pc_panel_plot"
  ),
  function_path
)


#----
#---- Load ERCCIS SQLite DB and tbl
#----
sqlite_file <- paste0(shareDir, "erccis_oct2020.sqlite")
table_name <- "erccis_oct2020"
erccis_data <- erccis.db.tbl(
  sqlite_file,
  table_name,
  sppGrps = sppNameTab$SpeciesGroup,
  sppNameLookUp = sppNameTab
)

list.files('Outputs/Temporal_coverage_R1_continuous/', pattern) 

#----
#---- Run cover estimates for continuous accumulative time and 
#---- decadal moving window
#----

#-- Continuous
run_cover(
  sppLookUp = sppNameTab,
  yearStart = 1960,
  grains = 1,
  yearType = "continuous",
  outPath = 'Outputs/Temporal_coverage_R1_continuous/',
  collateOnly = FALSE
)

#-- Decadal
run_cover(
  sppLookUp = sppNameTab,
  yearStart = 1970,
  yearType = "decadal",
  outPath = 'Outputs/Temporal_coverage_R1_decadal/',
  collateOnly = FALSE
) 

#---- Plot coverage index and extract statistics


#-- Q0 continuous
ciContQ0 <-
  ci_panel_plot(
    filePath= 'Outputs/Temporal_coverage_R1_continuous/',
    outPath='Outputs/Figures_R1/CI_Cont_Q0',
    q = 'SC_q0'
  )

#-- Q1 continuous
ciContQ1 <-
  ci_panel_plot(
    'Outputs/Temporal_coverage_R1_continuous/',
    'Outputs/Figures_R1/CI_Cont_Q1',
    q = 'SC_q1'
  )

#-- Q0 decadal
ciDecadalQ0  <-
  ci_panel_plot(
    'Outputs/Temporal_coverage_R1_decadal/',
    'Outputs/Figures_R1/CI_Decadal_Q0',
    q = 'SC_q0'
  )


#-- Q1 decadal
ciDecadalQ1 <-
  ci_panel_plot(
    'Outputs/Temporal_coverage_R1_decadal/',
    'Outputs/Figures_R1/CI_Decadal_Q1',
    q = 'SC_q1'
  )

#---- Plot proportion sites reaching completeness (i.e. >80%) and extract 
#---- statistics


#-- Q0 continuous
pc80ContQ0 <-
  pc_panel_plot(
    'Outputs/Temporal_coverage_R1_continuous/',
    'Outputs/Figures_R1/PC80_Cont_Q0',
    threshold = 0.8,
    q = 'SC_q0'
  )

#-- Q1 continuous
pc80ContQ1 <-
  pc_panel_plot(
    'Outputs/Temporal_coverage_R1_continuous/',
    'Outputs/Figures_R1/PC80_Cont_Q1',
    threshold = 0.8,
    q = 'SC_q1'
  )

#-- Q0 decadal
pc80DecadalQ0  <-
  pc_panel_plot(
    'Outputs/Temporal_coverage_R1_decadal/',
    'Outputs/Figures_R1/PC80_Decadal_Q0',
    threshold = 0.8,
    q = 'SC_q0'
  )

#-- Q1 decadal
pc80DecadalQ1 <-
  pc_panel_plot(
    'Outputs/Temporal_coverage_R1_decadal/',
    'Outputs/Figures_R1/PC80_Decadal_Q1',
    threshold = 0.8,
    q = 'SC_q1'
  )

#----
#---- Map 2019 sample completeness agreement 
#----

#-- Extract site / grain specific sample coverage estimates
sc_cont_Q1 <-
  propCovThres2('Outputs/Temporal_coverage_R1_continuous/',
                q = 'SC_q1',
                max_pc = TRUE) %>%
  replace_na(list(estCov  = 0)) %>%
  group_by(sppGrp, site, grain) %>%
  summarise(covMax = max(estCov, na.rm = TRUE)) 
sc_deca_Q1 <-
  propCovThres2('Outputs/Temporal_coverage_R1_decadal/',
                q = 'SC_q1',
                max_pc = TRUE) %>%
  replace_na(list(estCov  = 0)) %>%
  filter(year == 2019) %>%
  group_by(sppGrp, site, grain) %>%
  summarise(covMax = max(estCov, na.rm = TRUE))

#-- Count number of groups per site exceeding 80% coverage
n_groups_cont_Q1 <- sc_cont_Q1 %>%
  mutate(covComp = ifelse(covMax >= 0.8, 1, 0)) %>%
  group_by(site, grain) %>%
  summarise(nSppGrp = sum(covComp)) %>%
  rename(id = site) %>%
  mutate(id = as.character(id), time = "Continuous")
n_groups_deca_Q1 <- sc_deca_Q1 %>%
  mutate(covComp = ifelse(covMax >= 0.8, 1, 0)) %>%
  group_by(site, grain) %>%
  summarise(nSppGrp = sum(covComp)) %>%
  rename(id = site) %>%
  mutate(id = as.character(id), time = "Decadal")
n_groups_Q1 <- rbind(n_groups_cont_Q1,
                     n_groups_deca_Q1)

#-- Join counts to hex grids
hex1 <- left_join(hexGrain1, n_groups_Q1, by = c("id", "grain"))
hex4 <- left_join(hexGrain4, n_groups_Q1, by = c("id", "grain"))
hex16 <- left_join(hexGrain16, n_groups_Q1, by = c("id", "grain"))
hex64 <- left_join(hexGrain64, n_groups_Q1, by = c("id", "grain"))
covComCon <- do.call(rbind, list(hex1, hex4, hex16, hex64))
save(covComCon,
     file = "Outputs/Figures_R1/Coverage_completeness_n_species_Q1.Rdata")
#load("Outputs/Coverage_completeness_n_species.Rdata")

#-- Summarize
sumCovCom <- covComCon %>%
  st_drop_geometry() %>%
  group_by(grain, time) %>%
  summarise(medianGrp = median(nSppGrp), maxGrp = max(nSppGrp))
covComCon %>%
  st_drop_geometry() %>%
  group_by(grain, time) %>%
  mutate(N = n()) %>%
  group_by(grain, nSppGrp, time, N) %>%
  summarise(n_cells = n()) %>%
  mutate(pc = n_cells / N) %>%
  group_by(grain, time) %>%
  arrange(grain, desc(nSppGrp)) %>%
  mutate(pc_cum = cumsum(pc), n_cell_cum = cumsum(n_cells)) %>%
  print(n=100)

#-- Set the factor levels
covComCon$nSppGrp <- factor(covComCon$nSppGrp, levels = 0:15)
covComCon$time[covComCon$time == "Continuous"] <- "1960-2019"
covComCon$time[covComCon$time == "Decadal"] <- "2010-2019"
covComCon$time <-
  factor(covComCon$time, levels = c("1960-2019", "2010-2019"))
#-- Names for labelling facets
covComCon$grain[covComCon$grain == 1] <- '"Grain: 1 km"^2'
covComCon$grain[covComCon$grain == 4] <- '"Grain: 4 km"^2'
covComCon$grain[covComCon$grain == 16] <- '"Grain: 16 km"^2'
covComCon$grain[covComCon$grain == 64] <- '"Grain: 64 km"^2'
covComCon$grain <- factor(
  covComCon$grain,
  levels = c(
    '"Grain: 1 km"^2' ,
    '"Grain: 4 km"^2',
    '"Grain: 16 km"^2',
    '"Grain: 64 km"^2'
  )
)

#-- Small UK insert map
uk_counties <-
  st_read(
    paste0(
      shareDir,
      "/BoundaryLineUK/boundaryline_3962876/Boundary-line-historic-counties_region.shp"
    )
  )[, 1]
uk_counties$f <- 0
uk_counties$f[uk_counties$NAME == "Cornwall"] <- 1
uk_counties$f <- factor(uk_counties$f, levels = c(0, 1))
uk_map <- ggplot(uk_counties[, 3]) +
  geom_sf(aes(fill = f)) +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  theme_map() %+replace% theme(legend.position = 'none')

#-- Create plot
covComCon_Panel <- ggplot() +
  geom_sf(
    aes(),
    fill = 'transparent',
    colour = 'Black',
    size = 0.5,
    alpha = 0.5,
    data = cornwall_sf
  ) +
  geom_sf(
    aes(fill = nSppGrp),
    colour = 'transparent',
    size = 0.01,
    alpha = 0.75,
    data = covComCon
  ) +
  scale_colour_viridis(name = 'Number\nspp grps', discrete = TRUE) +
  scale_fill_viridis(name = 'Number\nspp grps', discrete = TRUE) +
  theme_bw(base_size = 14) %+replace%
  theme(legend.position = 'right',
        strip.background = element_rect(fill = "white")) +
  facet_grid(time ~ grain,  labeller = label_parsed)

#-- Add UK insert with Cornwall highlighted
corn_coord <- as.data.frame(st_coordinates(cornwall_sf))
x_min <- min(corn_coord$X)
x_max <- (max(corn_coord$X) + x_min) / 2
y_max <- max(corn_coord$Y)
y_min <- (min(corn_coord$Y) + y_max) / 2
annotation_custom2 <-
  function (grob,
            xmin = -Inf,
            xmax = Inf,
            ymin = -Inf,
            ymax = Inf,
            data) {
    layer(
      data = data,
      stat = StatIdentity,
      position = PositionIdentity,
      geom = ggplot2:::GeomCustomAnn,
      inherit.aes = TRUE,
      params = list(
        grob = grob,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      )
    )
  }
covComCon_Panel <- covComCon_Panel +
  annotation_custom2(
    ggplotGrob(uk_map),
    xmin = 110000,
    xmax = x_max,
    ymin = 45000,
    ymax = y_max,
    data = data.frame(grain = '"Grain: 64 km"^2', time = "2010-2019")
  )
ggsave(
  filename = "Outputs/Figures_R1/Spp_Sampling_Coverage_80_Panel_Q1.png",
  plot = covComCon_Panel,
  width = 15,
  height = 8
)

#----
#---- Correlation between coverage and sample number
#---- 

#-- Calculate correlation for all groups/grains
scMax_cont_Q1 <-
  maxEstCov(filePath = 'Outputs/Temporal_coverage_R1_continuous/',
            q = "SC_q1")
scMax_cont_Q1 %>% filter(grain == 1, covMax > 0)

covRecCorAll <-
  lapply(as.character(unique(sppNameTab[, 2])), function(x) {
    res <- lapply(c(1, 4, 16, 64), function(j) {
      tmp <- scMax_cont_Q1 %>% filter(grain == j & sppGrp == x)
      sppGrnCor <- cor(tmp$covMax, tmp$nRec, method = c("pearson"))
      resSpp <-
        data.frame('sppGrp' = x,
                   'grain' = j,
                   'corCovRec' = sppGrnCor)
    })
    resOut <- do.call(bind_rows, res)
  })
covRecCorAll <- do.call(bind_rows, covRecCorAll)
covRecCorAll %>%
  group_by(grain) %>%
  summarise(
    corGrn = mean(corCovRec),
    corGrnMin = min(corCovRec),
    corGrnMax = max(corCovRec)
  )
corAll <-
  covRecCorAll %>% pivot_wider(id_cols = 'sppGrp',
                               names_from = 'grain',
                               values_from = 'corCovRec')
write.csv(corAll,
          "Outputs/Figures_R1/Correlation_coeff_covXnRec_Q1.csv",
          row.names = FALSE)

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

#-------------------------------------------------------------------------------
#
# SECTION 5. Process environmental data for later analysis (clip, project, join)
#
#-------------------------------------------------------------------------------

#----
#---- Roads and pathways
#----
hwDatPath <-
  paste0(
    shareDir,
    'HighwaysAll_1641946/MasterMap Highways Network_3818259/Highways_Network_March19.gdb/Highways_Network_March19.gdb'
  )
hwDatLayers <- sf::st_layers(dsn = hwDatPath)

#--- Pathways
pathDat <- sf::st_read(dsn = hwDatPath, layer = "PathLink_N")
#- Need to shrink dataset by cropping cornwall
#- First need to transform the cornwall boundary to same projection
#- as the pathway data
cw <- cornwall_sf %>% st_transform(crs = st_crs(pathDat))
#- Then, crop the pathway data by the cornwall boundary
pathDatCw <- pathDat %>% st_crop(cw)
st_write(pathDatCw, paste0(shareDir, 'Cornwall_Pathways.shp'))
pathCw <- st_read(paste0(shareDir, 'Cornwall_Pathways.shp')) %>%
  select(TOID, geometry)
plot(pathCw[, 1])

#--- Roads
roadDat <- sf::st_read(dsn = hwDatPath, layer = "RoadLink_N")
cw <- cornwall_sf %>% st_transform(crs = st_crs(roadDat))
roadDatCw <- roadDat %>% st_crop(cw)
st_write(roadDatCw, paste0(shareDir, 'Cornwall_Roads.shp'))
roadCw <- st_read(paste0(shareDir, 'Cornwall_Roads.shp')) %>%
  select(TOID, geometry)

#--- Roads
linkDat <- sf::st_read(dsn = hwDatPath, layer = "ConnectingLink_N")
cw <- cornwall_sf %>% st_transform(crs = st_crs(linkDat))
linkDatCw <- linkDat %>% st_crop(cw)
st_write(linkDatCw, paste0(shareDir, 'Cornwall_Links.shp'))
linkCw <- st_read(paste0(shareDir, 'Cornwall_Links.shp'))  %>%
  select(TOID, geometry)

#--- Combine all
hwAll <- do.call(bind_rows, list(pathCw, roadCw, linkCw))
hwAllPrj <- hwAll %>% st_transform(crs = st_crs(cornwall_sf))
st_write(hwAllPrj, paste0(shareDir, 'Cornwall_Highways_All.shp'))
plot(hwAll[, 1])

#----
#---- Green spaces and protected areas
#----

#- Including greenspace, SSSI, National Parks, Nature reserves

#--- SSSI
SSSI <-
  st_read(
    paste0(
      shareDir,
      'NE_SSSI_Full/data/Sites_of_Special_Scientific_Interest_England.shp'
    )
  )
SSSI <- SSSI %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_crop(cornwall_sf)  %>%
  st_intersection(cornwall_sf) %>%
  select('geometry')
plot(SSSI[, 1])

#---  Greenspace
grnSpc <-
  st_read(
    paste0(
      shareDir,
      'GrnSpc2020_1641023/open-greenspace_3816376/GB_GreenspaceSite.shp'
    )
  )
grnSpc <- grnSpc %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_crop(cornwall_sf) %>%
  st_intersection(cornwall_sf) %>%
  select('geometry')
plot(grnSpc[, 1])

#--- National Nature Reserves
NNR <-
  st_read(paste0(
    shareDir,
    'NNR_Eng/National_Nature_Reserves___Natural_England.shp'
  ))
NNR <- NNR %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_crop(cornwall_sf) %>%
  st_intersection(cornwall_sf)
select('geometry')
plot(NNR[, 1])

#--- Special Areas of Conservation
SAC <-
  st_read(
    paste0(
      shareDir,
      'SAC_UK/Special_Areas_of_Conservation__England____Natural_England.shp'
    )
  )
SAC <- SAC %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_crop(cornwall_sf) %>%
  st_intersection(cornwall_sf) %>%
  select('geometry')
plot(SAC[, 1])

#--- Ancient woodland
AncWdLd <-
  st_read(
    paste0(
      shareDir,
      'Ancient_Woodland__England/Ancient_Woodland__England____Natural_England.shp'
    )
  )
AncWdLd <- AncWdLd %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_crop(cornwall_sf) %>%
  st_intersection(cornwall_sf) %>%
  select('geometry')
plot(AncWdLd[, 1])

#--- Combine all these polygons (if overlaps then won't be counted
#--- more than once)
paLd <-
  st_make_valid(do.call(bind_rows, list(SSSI, grnSpc, NNR, SAC, AncWdLd)))
paLd_union <- st_union(paLd) #%>% st_intersection(cornwall_sf)
st_write(paLd_union,
         paste0(shareDir, 'Cornwall_GrnSpc_SSSI_NNR_SAC_AncWdLd_All.shp'))
plot(paLd_union)

#~# #~# #~# Distance decay from major urban centres
urbanDD <-
  st_as_sf(
    readOGR(
      paste0(
        shareDir,
        '/strtgi_essh_gb/strtgi_essh_gb/data/urban_region.shp'
      )
    )
  ) %>%
  filter(CODE == 5420) %>%
  st_transform(crs = st_crs(cornwall_sf))
st_write(urbanDD , paste0(shareDir, 'Cornwall_Urban_Centres.shp'))
urbanArea_16 <- distance.decay(urbanArea, cornwall_sf, 16)

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

#-------------------------------------------------------------------------------
#
# SECTION 6. Process environmental predictors at appropriate grains
#
#-------------------------------------------------------------------------------

source.keep(c("feature.length","distance.decay",
              "feature.cover","ceh.landcov.extract"), function_path)

#----
#---- Extract at each grain
#----

#-- Highway data
highway <- st_read(paste0(shareDir, 'Cornwall_Highways_All.shp')) %>%
  st_transform(crs = st_crs(cornwall_sf))
#-- Urban data
urban <- st_read(paste0(shareDir, 'Cornwall_Urban_Centres.shp')) %>%
  st_transform(crs = st_crs(cornwall_sf)) 
#-- Green space data
grnSpc <-
  st_read(paste0(shareDir, 'Cornwall_GrnSpc_SSSI_NNR_SAC_AncWdLd_All.shp')) %>%
  st_transform(crs = st_crs(cornwall_sf)) %>%
  st_make_valid()

#-- 64 km2
highway_64 <-
  feature.length(highway, cornwall_sf, 64, hexGrain64, 
                 var_name = 'highway')
urbanDD_64 <-
  distance.decay(urban, cornwall_sf, 64, hexGrain64, 
                 var_name = 'urbanDD')
grnSpc_64 <-
  feature.cover(grnSpc, cornwall_sf, 64, hexGrain64, 
                var_name = 'grnSpc')
enviro_64 <-
  do.call(bind_rows, list(highway_64, urbanDD_64, grnSpc_64))

#-- 16 km2
highway_16 <-
  feature.length(highway, cornwall_sf, 16, hexGrain16,
                 var_name = 'highway')
urbanDD_16 <-
  distance.decay(urban, cornwall_sf, 16, hexGrain16,
                 var_name = 'urbanDD')
grnSpc_16 <-
  feature.cover(grnSpc, cornwall_sf, 16, hexGrain16,
                var_name = 'grnSpc')
enviro_16 <-
  do.call(bind_rows, list(highway_16, urbanDD_16, grnSpc_16))

#-- 4 km2
highway_4 <-
  feature.length(highway, cornwall_sf, 4, hexGrain4, 
                 var_name = 'highway')
urbanDD_4 <-
  distance.decay(urban, cornwall_sf, 4, hexGrain4, 
                 var_name = 'urbanDD')
grnSpc_4 <-
  feature.cover(grnSpc, cornwall_sf, 4, hexGrain4, 
                var_name = 'grnSpc')
enviro_4 <- do.call(bind_rows, list(highway_4, urbanDD_4, grnSpc_4))

#-- 1 km2
highway_1 <-
  feature.length(highway, cornwall_sf, 1, hexGrain1, 
                 var_name = 'highway')
urbanDD_1 <-
  distance.decay(urban, cornwall_sf, 1, hexGrain1, 
                 var_name = 'urbanDD')
grnSpc_1 <-
  feature.cover(grnSpc, cornwall_sf, 1, hexGrain1, 
                var_name = 'grnSpc')
enviro_1 <- do.call(bind_rows, list(highway_1, urbanDD_1, grnSpc_1))

# enviro_16 <- enviro_all %>% filter(grain == 16)
# enviro_64 <- enviro_all %>% filter(grain == 64)
# enviro_256 <- enviro_all %>% filter(grain == 256)
# tst <- left_join(hex1, highway_1, by = "id")
# ggplot(tst) + geom_sf(aes(fill=value)) + scale_fill_viridis()


#~# Environment at all grains
enviro_all <-
  do.call(bind_rows,
          list(enviro_1, enviro_4, enviro_16, enviro_64))
enviro_all$value[is.na(enviro_all$value)] <-
  0 # I forgot to set NAs to zeros earlier
save(enviro_all,
     file = 'Outputs/Environmental_variables_grain_R1.Rdata',
     compress = 'xz')

#~# Tidy workplace
rm(list = ls()[!ls() %in% c(keep_files)])

#-------------------------------------------------------------------------------
#
# SECTION 7. Sampling statistics (N observers, Scheme diversity)
#
#-------------------------------------------------------------------------------

#' Here, we extract statistics per grid cell on the:
#' 1) Number of unique recorders (mean across years)
#' 2) Some measure of scheme diversity

source.keep(c("erccis.db.tbl","species.group.select",
              "simpsonIndex","recToHexSetUp", "db.transform.proj"), function_path)

#---- Load ERCCIS SQLite DB and tbl
sppNameTab$SpeciesGroup[7] <- "Terrestrial mammal - bat"
erccis_data <- erccis.db.tbl(
  sqlite_file = paste0(shareDir, "erccis_structure_Nov2021.sqlite"),
  table_name = "erccis_structure_Nov2021",
  sppGrps = sppNameTab$SpeciesGroup,
  sppNameLookUp = sppNameTab
)

#---- Scheme structural diversity by species group 1960-2019
for (sppGrpN in unique(sppNameTab[, 2])[7]) {
  print(sppGrpN)
  
  #- Subset to just species group
  spp_grp <-
    as.character(sppNameTab[, 1][which(sppNameTab[, 2] == sppGrpN)])
  erccis_sppGrp <- species.group.select(erccis_data,
                                        spp_grp,
                                        cornwall_sf,
                                        "Cornwall",
                                        sppSelect = NULL)
  
  #- Extract records per grain grn = 64
  for (grn in c(1, 4, 16, 64)) {
    hexGrain <- get(paste0("hexGrain", grn))
    sppHex <-
      recToHexSetUp(
        bioDat = erccis_sppGrp,
        timeSlice = 1,
        region_sf = cornwall_sf,
        hexCellArea = grn,
        hexGrid = hexGrain
      )
    sppHex <- rename(sppHex,
                     sppGrp = SpeciesGroupFile,
                     StSurv = Structure_Survey)
    
    #- Scheme diversity
    sppHex$StSurv[grepl(sppHex$StSurv,
                        pattern = "Ad Hoc")] <- "None"
    sppHex$StSurv[is.na(sppHex$StSurv)] <- "None"
    StSurvN <- sppHex %>%
      group_by(id, sppGrp, grain, StSurv) %>%
      dplyr::summarise(StSurv_N = n())
    schemeDiv <- StSurvN %>%
      group_by(id, sppGrp, grain) %>%
      summarise(D = simpsonIndex(StSurv_N)) %>%
      dplyr::select(id, sppGrp, grain, D) %>%
      st_drop_geometry()
    schemeN <- StSurvN %>%
      group_by(id, sppGrp, grain) %>%
      filter(StSurv != "None") %>%
      summarise(nSchemes = length(StSurv_N)) %>%
      dplyr::select(id, sppGrp, grain, nSchemes) %>%
      st_drop_geometry()
    
    #- Proportion of ad hoc (conservative)
    propAdhoc <- sppHex %>%
      group_by(id, sppGrp, grain, StSurv) %>%
      dplyr::summarise(StSurv_N = n()) %>%
      mutate(StSurv = ifelse(StSurv == "None", "US", "ST")) %>%
      group_by(id, sppGrp, grain) %>%
      mutate(N = sum(StSurv_N)) %>%
      filter(StSurv == "US") %>%
      mutate(prop_adhoc = StSurv_N / N) %>%
      dplyr::select(id, sppGrp, grain, prop_adhoc) %>%
      st_drop_geometry()
    
    #- Number of records
    nRecorders <- sppHex %>%
      group_by(id, sppGrp, grain, RecorderAlias) %>%
      dplyr::summarise(nRecorders = n()) %>%
      group_by(id, sppGrp, grain) %>%
      summarise(nRecorders = n()) %>%
      st_drop_geometry()
    
    #- Join
    scheme_metrics <-
      Reduce(function(...)
        left_join(...,
                  by = c("id", "sppGrp", "grain")),
        list(schemeDiv, schemeN, propAdhoc, nRecorders))
    
    #- Save
    save(
      scheme_metrics,
      file = paste0(
        "Outputs/Recording_scheme_metrics/SchemeMetrics_Spp_",
        sppGrpN,
        "_grain_",
        grn,
        ".Rdata"
      )
    )
    
  }
  
}


#---- Scheme information
f <-
  list.files("Outputs/Recording_scheme_metrics/", full.names = TRUE)
scheme_metrics <- lapply(f, function(x) {
  tmp <- get(load(x))
  tmp$sppGrp <- strsplit(x, split = "Spp_|_grain")[[1]][2]
  tmp$id <- as.numeric(tmp$id)
  tmp
})
scheme_metrics <- do.call(rbind, scheme_metrics) 
save(scheme_metrics,
     file = 'Outputs/scheme_metrics.Rdata',
     compress = 'xz')

#-- Number of schemes per cell
scheme_metrics %>%
  group_by(grain) %>%
  dplyr::summarise(N = mean(nSchemes, na.rm = TRUE),
                   l95  = quantile(nSchemes, prob = 0.025, na.rm = TRUE),
                   u95  = quantile(nSchemes, prob = 0.975, na.rm = TRUE))
  
#-- Proportion of opportunistic
scheme_metrics %>%
  group_by(grain) %>%
  dplyr::summarise(N = mean(prop_adhoc, na.rm = TRUE),
                   l95  = quantile(prop_adhoc, prob = 0.025, na.rm = TRUE),
                   u95  = quantile(prop_adhoc, prob = 0.975, na.rm = TRUE))

#-- Number of recorders
scheme_metrics %>%
  group_by(grain) %>%
  dplyr::summarise(N = mean(nRecorders, na.rm = TRUE),
                   Min = min(nRecorders, na.rm = TRUE),
                   Max = max(nRecorders, na.rm = TRUE))

#-- Diversity of records
scheme_metrics %>%
  group_by(grain) %>%
  dplyr::summarise(N = mean(D, na.rm = TRUE),
                   l95  = quantile(D, prob = 0.025, na.rm = TRUE),
                   u95  = quantile(D, prob = 0.975, na.rm = TRUE))




#-------------------------------------------------------------------------------
#
# SECTION 8. Spatial analysis of sampling completeness vs environmental vars
#
#-------------------------------------------------------------------------------

source.keep(c(
  "maxEstCov",
  "makeHexGrid",
  "prepareModelDat",
  "create.nb.mat",
  "modCheck",
  "bayesModsRun",
  "species.level.models",
  "sppModResExtr"
),
function_path)

#-- Environmental data
enviro_all <-
  get(load('Outputs/Environmental_variables_grain_R1.Rdata')) %>%
  pivot_wider(
    id_cols = c(id, grain),
    names_from = var_name,
    values_from = value
  ) %>%
  mutate(id = as.integer(id)) %>%
  group_by(grain) %>%
  rename(remoteness = urbanDD, 
         accessibility = highway,
         natureSpace = grnSpc) %>%
 mutate(across(c(accessibility, remoteness, natureSpace), scale)) 
  

#-- Scheme information
scheme_metrics <- get(load('Outputs/scheme_metrics.Rdata')) %>%
  mutate(nRecLog10 = log10(nRecorders)) %>% 
  group_by(sppGrp, grain) %>%
  mutate(across(nRecLog10, scale))

#-- Load max coverage per gain for all species groups
maxCov <- maxEstCov('Outputs/Temporal_coverage_R1_continuous/')
maxCov <- maxCov %>% rename('id' = 'site')

#-- Join maxCov with environmental data
maxEnv <- left_join(maxCov, enviro_all, by = c('id', 'grain'))
#-- Join scheme metrics
maxEnv <- left_join(maxEnv, scheme_metrics, by = c('id', 'grain', 'sppGrp'))

#---
#--- Grain 1 
#---

#- Model data
DAT1 <- prepareModelDat(hexGrain1, 1, maxEnv)
#- Create neighbourhood matrix 3 deep
W <- create.nb.mat(hexGrain1, n = 24, DAT1) # this is all cells within three

#- Fit model
mod1km <-
  brm(
    covComp |
      weights(area) ~ remoteness + accessibility + 
      natureSpace + nRecLog10 + (1|sppGrp),# + car(W, id, type = 'escar'),
    data = DAT1,
    data2 = list(W = W[[2]]),
    family = bernoulli(),
    chain = 2,
    iter = 12000,
    thin = 5,
    cores = 2,
    prior=set_prior(horseshoe(df = 3, par_ratio = 0.1)),
    control = list(adapt_delta = 0.95)
  )
save(mod1km, file = "Outputs/Figures_R1/brms_models/brms_CAR_1km.Rdata",
     compress = 'xz')
unique(mod1km$data$id)
#- Check for spatial autocorrelation and model fit
modCheck(mod1km, DAT1)
summary(mod1km)
bayes_R2(mod1km)
pp_check(mod1km, resp = "covComp")
#- Test specific hypothesis
hypothesis(mod1km, 'remoteness < 0')
hypothesis(mod1km, 'accessibility > 0')
hypothesis(mod1km, 'natureSpace > 0')
hypothesis(mod1km, 'nRecLog10 > 0')
nrow(DAT1)

#---
#--- Grain 4
#---

#- Model data
DAT4 <- prepareModelDat(hexGrain4, 4, maxEnv)
#- Create neighbourhood matrix 3 deep
W <-
  create.nb.mat(hexGrain4, n = 12, DAT4) # this is all cells within three

#- Fit model
mod4km <-
  brm(
    covComp |
      weights(area) ~ remoteness + accessibility + 
      natureSpace + nRecLog10 + (1|sppGrp) + 
      car(W, id, type = 'escar'),
    data = DAT4,
    data2 = list(W = W[[2]]),
    family = bernoulli(),
    chain = 2,
    iter = 12000,
    thin = 10,
    cores = 2,
    control = list(adapt_delta = 0.95)
  )
save(mod4km, 
     file = "Outputs/Figures_R1/brms_models/brms_car_4km.Rdata",
     compress = 'xz')
#- Check for spatial autocorrelation and model fit
modCheck(mod4km, DAT4)
summary(mod4km)
bayes_R2(mod4km)
pp_check(mod4km, resp = "covComp")
#- Test specific hypothesis
#- Test specific hypothesis
hypothesis(mod4km, 'remoteness < 0')
hypothesis(mod4km, 'accessibility > 0')
hypothesis(mod4km, 'natureSpace > 0')
hypothesis(mod4km, 'nRecLog10 > 0')
nrow(DAT4)

#---
#--- Grain 16 
#---

#- Model data
DAT16 <- prepareModelDat(hexGrain16, 16, maxEnv)
#- Create neighbourhood matrix 3 deep
W <- create.nb.mat(hexGrain16, n = 8, DAT16) # this is all cells within two

#- Fit model
mod16km <-
  brm(
    covComp |
      weights(area) ~ remoteness + accessibility + 
      natureSpace + nRecLog10 + (1|sppGrp) + 
      car(W, id, type = 'escar'),
    data = DAT16,
    data2 = list(W = W[[2]]),
    family = bernoulli(),
    chain = 2,
    iter = 12000,
    thin = 10,
    cores = 2,
    prior=set_prior(horseshoe(df = 3, par_ratio = 0.1)),
    control = list(adapt_delta = 0.95)
  )
save(
  mod16km,
  file = "Outputs/Figures_R1/brms_models/brms_car_16km2.Rdata",
  compress = 'xz'
)
load("Outputs/Figures_R1/brms_models/brms_car_16km.Rdata")
#- Check for spatial autocorrelation and model fit
modCheck(mod16km, DAT16)
summary(mod16km)
bayes_R2(mod16km)
pp_check(mod16km, resp = "covComp")
#- Test specific hypothesis
hypothesis(mod16km, 'remoteness < 0')
hypothesis(mod16km, 'accessibility > 0')
hypothesis(mod16km, 'natureSpace > 0')
hypothesis(mod16km, 'nRecLog10 > 0')


#~# #~# #~#
#~# #~# #~# Fit species level models #~# #~# #~#

#~# 1 km 
species.level.models(hexCellArea = 1, hexGrain1, covEnvDat = maxEnv, 15, unique(sppNameTab[, 2]))

#~# 4 km 
species.level.models(hexCellArea = 4, hexGrain4, covEnvDat = maxEnv, 12, unique(sppNameTab[, 2]))

#~# 16 km 
species.level.models(hexCellArea = 16, hexGrain16, covEnvDat = maxEnv, 8, unique(sppNameTab[, 2]))

#~# 64 km 
#species.level.models(region_sf = cornwall_sf, hexCellArea = 64, covEnvDat = maxEnv)

#~# #~# #~# Species level covariate plots #~# #~# #~#

#~# Extract results
sppModsRes <- sppModResExtr(sppGrp = unique(sppNameTab[,2]), grains = c(16))
sppModsRes <- sppModsRes[,!grepl(names(sppModsRes), pattern="Intercept")]

# Median Estimate
dat_b <- sppModsRes  %>%
  select(sppGrp, grain, nObs, grep(names(sppModsRes), pattern="_b")) %>%
  pivot_longer(cols = c(ends_with("b")), names_to = 'covariate', values_to = 'Estimate') %>%
  filter(covariate != "R2_b")
dat_b$covariate <- sapply(dat_b$covariate, function(x) strsplit(x, split = "_")[[1]][1] )

# Lower CI
dat_L <- sppModsRes  %>%
  select(sppGrp, grain, nObs, grep(names(sppModsRes), pattern="_L")) %>%
  pivot_longer(cols = c(ends_with("L")), names_to = 'covariate', values_to = 'Lower') %>%
  filter(covariate != "R2_L")
dat_L$covariate <- sapply(dat_L$covariate, function(x) strsplit(x, split = "_")[[1]][1] )

# Upper CI
dat_U <- sppModsRes  %>%
  select(sppGrp, grain, nObs, grep(names(sppModsRes), pattern="_U")) %>%
  pivot_longer(cols = c(ends_with("U")), names_to = 'covariate', values_to = 'Upper') %>%
  filter(covariate != "R2_U")
dat_U$covariate <- sapply(dat_U$covariate, function(x) strsplit(x, split = "_")[[1]][1] )
  
# Hypothesis
dat_H1 <- sppModsRes  %>%
  select(sppGrp, grain, nObs, grep(names(sppModsRes), pattern="_H1_PProb")) %>%
  pivot_longer(cols = c(ends_with("H1_PProb")), names_to = 'covariate', values_to = 'H1') %>%
  filter(covariate != "R2_H1")
dat_H1$covariate <- sapply(dat_H1$covariate, function(x) strsplit(x, split = "_")[[1]][1] )

# Combine
dat_long <-
  Reduce(function(...)
    left_join(..., by = c('sppGrp', 'grain', 'nObs', 'covariate')),
    list(dat_b, dat_L, dat_U, dat_H1))
dat_long$covariate <-
  factor(
    dat_long$covariate,
    levels = c(
      "Remoteness",
      "Accessibility",
      "NatureSpace",
      "NumberOfRecorders"
    )
  )

# Evidence ratio
dat_long$H1PP <- '0'
dat_long$H1PP[dat_long$H1 > 0.9] <- '1'

# Facet labels
dat_long$labfacet <- NA
dat_long$labfacet[dat_long$grain == 1] <- '"Grain: 1 km"^2'
dat_long$labfacet[dat_long$grain == 4] <- '"Grain: 4 km"^2'
dat_long$labfacet[dat_long$grain == 16] <- '"Grain: 16 km"^2'
dat_long$labfacet <-
  factor(dat_long$labfacet,
         levels = c('"Grain: 1 km"^2', '"Grain: 4 km"^2', '"Grain: 16 km"^2'))

#~# Create plot for 1 km (as 4km is very similary)
spCovPlot1km <- ggplot(dat_long) +
  geom_point(aes(x = sppGrp, y = Estimate, colour = H1PP), size = 2) +
  geom_linerange(aes(
    x = sppGrp,
    ymin = Lower,
    ymax = Upper,
    colour = H1PP
  ), size = 0.75) +
  scale_colour_manual(values = c('0' = 'Black', '1' = '#E41A1C'),
                      drop = F) +
  geom_hline(yintercept = 0,
             linetype = 2,
             colour = 'Black') +
  theme_bw(base_size = 14) %+replace% theme(axis.title.y = element_blank(), legend.position = 'none') +
  scale_x_discrete(labels = sub('_', ' ', unique(sppNameTab[, 2])), drop = F) +
  ylab("Standardised regression coefficients") +
  coord_flip() + facet_grid(labfacet ~ covariate, labeller = label_parsed, scale = 'free')
ggsave(
  filename = 'Outputs/Figures_R1/Figure_4_Species_level_covariate_Grains_all.png',
  plot = spCovPlot1km,
  width = 12,
  height = 8
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SECTION 9. Use Zeta diversity to check for presence of turnover through time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source.keep(c(
  "maxEstCov",
  "erccis.db.tbl",
  "siteXSpeciesXTime",
  "species.group.select",
  "temporal.zeta",
  "db.transform.proj",
  "recToHexSetUp"
),
function_path)

#----
#---- Load ERCCIS SQLite DB and tbl
#----
sqlite_file <- paste0(shareDir, "erccis_oct2020.sqlite")
table_name <- "erccis_oct2020"
erccis_data <- erccis.db.tbl(
  sqlite_file,
  table_name,
  sppGrps = sppNameTab$SpeciesGroup,
  sppNameLookUp = sppNameTab
)


#-- Load max coverage per gain for all species groups
maxCov <- maxEstCov('Outputs/Temporal_coverage_R1_continuous/')
maxCov <- maxCov %>% rename('site' = 'id')

#~# #~# #~# Estimated coverage #~# #~# #~#
for (sppGrpN in unique(sppNameTab[, 2])[2:15]) {
  print(sppGrpN)
  #~# Subset to just species group
  erccis_sppGrp <-
    species.group.select(erccis_data, 
                         sppNameTab[, 1][which(sppNameTab[, 2] == sppGrpN)], 
                         cornwall_sf, 
                         sppSelect = NULL)
  
  #~# Grains
  for (grn in c(4,16,64)) {
    #~# Get pre-made hex grid
    hexGrain <- get(paste0('hexGrain', grn))
    
    #~# Cover estimate
    zeta_tmp <- temporal.zeta(
      erccis_sppGrp,
      sppGrpN,
      timeSlice = 5,
      region_sf = cornwall_sf,
      grain = grn,
      hexGrid = hexGrain
    )
    
    write.csv(
      zeta_tmp,
      file = paste0(
        'Outputs/Figures_R1/Temporal_zeta_decay/',
        sppGrpN,
        '_',
        grn,
        '.csv'
      ),
      row.names = FALSE
    )
    
  }
  
}
zeta_dd_all <-
  lapply(list.files('Outputs/Figures_R1/Temporal_zeta_decay/', full.names = TRUE), function(x)
    read.csv(x))
zeta_dd_all <- do.call(bind_rows, zeta_dd_all)
#~# Join maxCov with environmental data
maxEnv <-
  left_join(maxCov, zeta_dd_all, by = c('site', 'grain', 'sppGrp'))
#'
pc_dd_sim <- maxEnv %>% filter(!is.na(covMax)) %>%
  mutate(sig = if_else(dd_p < 0.01, 1, 0)) %>%
  group_by(sppGrp, grain) %>%
  summarise(n_sig = (sum(sig, na.rm = TRUE) / n()) * 100) %>% print(n = 60)

pc_dd_sim %>% filter(n_sig > 10)

pc_dd_sim %>% group_by(grain) %>% summarise(
  n_sig_m = median(n_sig),
  n_sig_min = min(n_sig),
  n_sig_max = max(n_sig)
)

#------------------------------------------------------------------------------#
#----------------------------------- END --------------------------------------#
#------------------------------------------------------------------------------#
