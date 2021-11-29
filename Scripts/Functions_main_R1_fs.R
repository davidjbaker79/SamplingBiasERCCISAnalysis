#-------------------------------------------------------------------------------
#
#     FUNCTIONS for Multidimensional analysis of local biodiversity records 
#      (using ERCCIS data as an example) - Revision 1
#       By David Baker (13/10/2021)
#
#-------------------------------------------------------------------------------

#----
#----  Load / Process ERCCISS data        
#----  

#' Load and, if necessary, install packages automatically
#'
#' Just saves a bunch of time, especially when moving between
#' computers.
#'
#' @param packageNames A vector of quoted package names.
#' @return A table of packages installed and loaded.
#' @export
package.load <- function(packageNames) {
  
  packageInstalled <- list.files(.libPaths())
  
  l <- lapply(packageNames, function(x) {
    
    if (!(x %in% packageInstalled)) utils::install.packages(x)
    
    require(x, character.only = TRUE)
    
  })
  
  report <- data.frame(Package = packageNames, LoadStatus = unlist(l))
  
}

#' Source only certain functions
#'
#' This helps keep on those functions needed for a particular
#' section of a long analysis load.
#'
#' @param keepFunctions Functions wanted in memory
#' @param functionPath Path to where functions are stored
#' @export
source.keep <- function(keepFunctions, functionPath) {
  
  # Load from source
  source(functionPath)
  
  # List functions
  loadedFunctions <- as.vector(utils::lsf.str(envir = globalenv()))
  
  # List functions to remove
  removeFunctions <- loadedFunctions[!(loadedFunctions %in% c(keepFunctions, "source.keep", "package.load"))]
  
  # Remove functions
  rm(list = removeFunctions, pos = ".GlobalEnv")
  
}

#' Load ERCCIS data
#'
#' Create a connection with the ERCCIS SQLite db and
#' filter missing coordinates and create required date formats
#'
#' @param sqlite_file A path name to .sqlite file.
#' @param table_name The name of .sqlite table to load (if multiple)
#' @param sppGrps Species group required.
#' @param sppNameLookUp A table with species names can be provided for 
#' filtering. A column can also
#' be added to add better names (e.g. for plotting).
#'
#' @import sf
#'
#' @return An tibble data.frame
#'
#' @export
erccis.db.tbl <- function(sqlite_file, table_name, 
                          sppGrps = NULL, sppNameLookUp = NULL) {
  
  # Create db connection
  sq_con <-
    RSQLite::dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_file)
  
  # Select species group
  qu <- paste0('SELECT * FROM ',table_name,' WHERE "SpeciesGroup" == :x')
  sq_dat <-
    RSQLite::dbGetQuery(sq_con,
                        qu,
                        params = list(x = sppGrps))
  if (is.null(sppGrps) | nrow(sq_dat) == 0) {
    qu <- paste0('SELECT DISTINCT SpeciesGroup FROM ', table_name)
    sppGrps <-
      RSQLite::dbGetQuery(sq_con, qu)
    cat('Please enter one or more of the following species group names')
    print(sppGrps)
  }
  
  # Remove records with no spatial coordinates (typically a small percentage)
  names(sq_dat)[which(names(sq_dat) == "SpatialLong")] <- 'X'
  names(sq_dat)[which(names(sq_dat) == "SpatialLat")] <- 'Y'
  sq_dat <- sq_dat[!is.na("X") & !is.na("Y"),]
  sq_dat <- sq_dat[sq_dat$SpeciesTaxonLevel == "Species", ]
  sq_dat <- sq_dat[!is.na(sq_dat$SpeciesGroup), ]
  sq_dat$StartDate <- as.Date(sq_dat$StartDate, format = "%Y-%m-%d")
  sq_dat$EndDate <- as.Date(sq_dat$EndDate, format = "%Y-%m-%d")
  
  # Add in better names if required
  if (!is.null(sppNameLookUp))
    sq_dat <-
    merge(sq_dat, sppNameLookUp, by = "SpeciesGroup", all.x = TRUE)
  sq_dat
  
}

#' Re-project point data
#'
#' This function uses a reference spatial polygon to transform the point data 
#' to a new projection
#' and returns an sf spatial object.
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and SpatialLat given as \code{Y}.
#' @param crs_sf_ref An sf object with the desired crs.
#' @param region A region defined within  \code{crs_sf_ref}.
#'
#' @import sf
#'
#' @return An sf object.
#'
#' @export
db.transform.proj <- function(bioDat, crs_sf_ref, 
                              region = NULL, df = FALSE) {
  
  # Convert to sf
  dat_sf <- st_as_sf(
    x = bioDat,
    coords = c("X", "Y"),
    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  )
  
  # Spatial transform and mask out records off the coast and Isle of Scilly (also few in Devon)
  dat_sf <- st_transform(dat_sf, crs = st_crs(crs_sf_ref))
  dat_sf <- st_join(dat_sf, crs_sf_ref, join = st_intersects)
  
  if (!is.null(region)) dat_sf <- dat_sf[dat_sf$name %in% region, ]
  
  if(df) {
    dat_df <- st_drop_geometry(dat_sf)
    dat_df$lon <- st_coordinates(dat_sf)[,1]
    dat_df$lat <- st_coordinates(dat_sf)[,2]
    dat_sf <- dat_df
  }
  dat_sf
  
}

#' Obtain species group data
#'
#' Wrapper function for filtering by species groups and projecting point data
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and SpatialLat given as \code{Y}.
#' @param sppGrp Name of species group (unwise to run on all species!!!).
#' @param region_sf An sf object defining the spatial region.
#' @param regionSelect A region defined within  \code{crs_sf_ref}.
#' @param sppSelect A SpeciesVenacular name or vector of names if just a subset 
#' of species is wanted.
#'
#' @return A filtered and transformed data frame.
#'
#' @export
species.group.select <- function(bioDat, sppGrp, region_sf, 
                                 regionSelect = NULL, sppSelect = NULL) {
  
  # Subset to just species group
  bioDat <- bioDat[bioDat$SpeciesGroup %in% sppGrp,]
  
  # If we only want a subset of species from a group
  if (!is.null(sppSelect))
    bioDat <- bioDat[bioDat$SpeciesVenacular %in% sppSelect,]
  
  # Create spatial object from db subset and transform to BNG
  bioDat <- db.transform.proj(bioDat, region_sf, regionSelect)
  
}

#' Plots for evaluating species distributions
#'
#' @param spp_dat An export from \code{db.transform.proj}
#' @param region_sf An sf object defining the spatial region.
#' @param focal_sp The SpeciesVenacular name.
#' @param start_yr A numeric year defining the start of the focal period.
#' @param end_yr A numeric year defining the end of the focal period.
#' @param plotType One of Map, Hist, or listLength.
#'
#' @import data.table
#' @import sf
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @return A ggplot figure.
#'
#' @export
#spp_dat = grpPrj; focal_sp =  "Hazel Dormouse"; start_yr = 2000; end_yr = 2019; region_sf = cornwall_sf; plotType = 'ListLength'
species.evaluation <- function(spp_dat,
                               focal_sp,
                               start_yr = 2001,
                               end_yr = NULL,
                               region_sf,
                               plotType = c("Map", "Hist", "ListLength")) {
  # Subset time period
  if (!is.null(start_yr) &
      !is.null(end_yr))
    spp_dat <- spp_dat[spp_dat$Yr >= start_yr & spp_dat$Yr <= end_yr,]
  if (plotType !=  "Map")
    spp_dat <-  data.table(spp_dat)
  fspp <- spp_dat[spp_dat$SpeciesVenacular == focal_sp,]
  
  # Plot target species occurrence data
  if (plotType == "Map") {
    plotOut <-
      ggplot() +
      geom_sf(fill = "grey75", data = region_sf) +
      geom_sf(
        aes_string(colour = "Yr"),
        alpha = 0.5,
        size = 1,
        data = fspp
      ) +
      scico::scale_colour_scico(palette = "bamako",
                                direction = -1,
                                name = "Year") +
      theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank()) +
      ggsn::north(region_sf, "topleft", 0.3, 16) +
      ggtitle(focal_sp)
    
  }
  
  # Temporal trend in regional record
  if (plotType == "Hist") {
    # Convert to data.table
    fspp_n <- fspp[, .N, by = list(Yr)]
    
    plotOut <-
      ggplot(aes_string(x = "Yr", y = "N"), data = fspp_n) +
      geom_col(fill = "grey75") +
      theme_bw(base_size = 14) +
      xlab("Year") +
      ylab("Number of Records") +
      ggtitle(focal_sp)
    
  }
  
  # List length per year
  if (plotType == "ListLength") {
    sp_grp <- spp_dat[, .N, by = list(SpatialReference, StartDate, Yr)]
    ll_n <- spGrp[, .N, by = list(Yr, N)]
    names(ll_n) <- c('Yr', 'll', 'N')
    
    
    sp_foc <- fspp[, .N, by = list(SpatialReference, StartDate, Yr)]
    foc_n <- sp_foc[, .N, by = list(Yr, N)]
    names(foc_n) <- c('Yr', 'll', 'N')
    
    plotOut <-
      ggplot() +
      geom_point(
        aes_string(x = "ll", y = "log10(N)", group = "Yr"),
        colour = "black",
        data = ll_n
      ) +
      geom_line(
        aes_string(x = "ll", y = "log10(N)", group = "Yr"),
        colour = "grey90",
        data = ll_n
      ) +
      geom_point(
        aes_string(x = "ll", y = "log10(N)", group = "Yr"),
        colour = "red",
        data = foc_n
      ) +
      theme_bw(base_size = 14) +
      xlab("List length") + ylab("Number of lists (log10)") +
      ggtitle(focal_sp)
    
  }
  
  plotOut
  
}


#' Make hexagonal grid
#'
#' Analysis is conducted on hexagonal grids and this function creates the grid 
#' with a specified area.
#'
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell 
#' (in km2).
#'
#' @import sf
#'
#' @return An sf object with 'id' used as a unique cell id.
#'
#' @export
makeHexGrid <- function(region_sf, hexCellArea) {
  
  hexGrid <-
    st_make_grid(region_sf,
                 cellsize = (1074.57 * sqrt(hexCellArea)),
                 square = FALSE)
  hexGrid <- st_as_sf(hexGrid)
  hexGrid$id <- rownames(hexGrid)
  hexGrid$grain <- hexCellArea
  hexGrid <- st_intersection(hexGrid, region_sf)
  hexGrid$area <- round(st_area(hexGrid) / 10 ^ 6, 2)
  hexGrid <- cbind(hexGrid, st_coordinates(st_centroid(hexGrid)))
  
}

#' Assign to hexagonal grid by any time slice
#'
#' Function assigns records to grid cells based on spatial overlap.
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and 
#' SpatialLat given as \code{Y}.
#' @param timeSlice A numeric value specifying the length on each timeslice 
#' (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell 
#' (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#'
#' @import sf
#'
#' @return A tibble data frame with hexagonal cell id and X, Y coordinates.
#'
#' @export
recToHexSetUp <- function(bioDat,
                          timeSlice,
                          region_sf,
                          hexCellArea,
                          hexGrid = NULL) {
  # Create time slices based on intervals
  start_date <- seq.Date(as.Date("1960/01/01"),
                         as.Date("2020/01/01"), "year")[seq(1, 61, timeSlice)]
  start_date <- start_date[!is.na(start_date)]
  
  # Create TS labels
  TSLab <- sub("-01-01", "", start_date)
  TSLab <- paste(TSLab[-length(start_date)], TSLab[-1], sep = "_")
  
  # For each time slice filter bioDat into time slices using start and end dates
  bioDat[["TS"]] <-
    cut(
      bioDat$StartDate,
      breaks = start_date,
      right = TRUE,
      labels = TSLab
    )
  bioDat <- bioDat[!is.na(bioDat$TS),]
  bioDat <-
    bioDat[, !(names(bioDat) %in% c("label", "name", "code"))]
  
  # Create hexagonal grid with id attribute
  if (is.null(hexGrid))
    hexGrid <- makeHexGrid(region_sf, hexCellArea)
  
  # Join biological records to hexGrid
  bioHex <- st_join(hexGrid, bioDat, join = st_intersects)
  
  # Add in XY
  idXY <-
    data.frame(st_coordinates(st_centroid(hexGrid)), id = hexGrid$id)
  bioHex <- merge(bioHex, idXY, by = c('id','X','Y'))
  
}

#' Number of records \code{x} hexagonal grid \code{x} time slice
#'
#' This function summarises number of records per cell and time slice.
#' This is useful for exploring trends in record collection for identifying 
#' periods when
#' records may be more or less representative.
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given 
#' as Y.
#' @param sppGrp Name of species group (unwise to run on all species!!!).
#' @param timeSlice A numeric value specifying the length on each timeslice 
#' (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell 
#' (in km2).
#' @param bioHex An export from \code{recToHexSetUp}.
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#'
#' @import sf
#' @import data.table
#' @importFrom rlang .data
#'
#' @return A data frame with hexagonal cell id and X, Y coordinates.
#'
#' @export
recXgridXtime <-
  function(bioDat,
           sppGrp,
           timeSlice,
           region_sf,
           hexCellArea,
           bioHex = NULL,
           hexGrid = NULL) {
    if (is.null(bioHex))
      bioHex <-
        recToHexSetUp(bioDat, timeSlice, region_sf, hexCellArea, hexGrid)
    
    # Summarise by id and TS
    nRecHex <- st_drop_geometry(bioHex)
    nRecHex <- data.table(nRecHex)
    nRecHex <- nRecHex[, .N, by = list(id, X, Y, TS)]
    nRecHex$sppGrp <- sppGrp
    nRecHex <- as.data.frame(nRecHex)
    
    return(nRecHex)
    
  }

#' Site \code{x} Species \code{x} Time data frame
#'
#' This function creates a site \code{x} species \code{x} time matrix
#' (presence/absence), which can then be subset into specific time slices.
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given 
#' as Y.
#' @param timeSlice A numeric value specifying the length on each timeslice 
#' (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea Area of each hexagonal grid cell.
#' @param bioHex An export from \code{recToHexSetUp}.
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param randTS If TRUE the date of surveys are randomised within each cell.
#'
#' @import data.table
#' @import sf
#' @importFrom rlang .data
#'
#' @return A binary (0,1) presence absence data frame is returned
#'
#' @export
siteXSpeciesXTime <- function(bioDat, timeSlice, region_sf, hexCellArea,
                              bioHex = NULL, hexGrid = NULL) {
  
  # Call setup function to create TS and hexGrid association
  if (is.null(bioHex))
    bioHex <-
      recToHexSetUp(bioDat, timeSlice, region_sf, hexCellArea, hexGrid)
  
  # Summarise
  sumHex <- st_drop_geometry(bioHex)
  sumHex <- data.table(sumHex)
  sumHex <-
    sumHex[, .N, by = list(id, TS, SpeciesScientific, SpatialReference)]
  sumHex$siteTime <- paste0(sumHex$id, "_", sumHex$TS, "_",
                            sumHex$SpatialReference)
  sumHex <- sumHex[, c("siteTime", "SpeciesScientific", "N")]
  sumHex <- as.data.frame(sumHex)
  
  # Convert to site x species df
  siteXSpp <- tidyr::pivot_wider(sumHex,
                                 names_from = SpeciesScientific,
                                 values_from = N)
  siteXSpp[is.na(siteXSpp)] <- 0
  siteXSpp <- as.data.frame(siteXSpp)
  
}

#~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~#
#~# #~# #~#         Sampling analysis               #~# #~# #~# 

#' Test sample size thresholds for coverage estimates
#'
#' @param exp_cond The combinations of experiments to run.
#' @param hexCellArea Area of each hexagonal grid cell.
#' @param hex1 An export from \code{recToHexSetUp}.
#' 
#' @import dplyr
#' @import fsSDM
#' @import vegan
#' @importFrom rlang .data
#'
#' @return A dataframe.
test_coverage_est <- function(exp_cond, hexGridArea, hex1) {
    # Create hex grids and associate 1km with coarser grid

    hex2 <- makeHexGrid(cornwall_sf, hexGridArea)
    hexRelate <- st_join(hex1, hex2) %>%
      st_drop_geometry() %>%
      dplyr::select(c(id.x, id.y)) %>%
      dplyr::rename(id = id.x, id2 = id.y)
    
    res_out <- lapply(1:nrow(exp_cond), function(i) {
      # Simulate multi-species occupancy
      spp_occ <-
        occupancy.multi(
          enviro_dat,
          enviro_vars,
          n_spp = 100,
          #exp_cond[i,1],
          pca_axes = 3,
          cornwall_sf,
          plotRichness = T,
          plotSpecies = F
        )
      
      # Generate virtual sampling
      surv_ex <-
        virtual.sampling(
          spp_occ,
          p_sites =  exp_cond[i, 2],
          n_yrs = 90,
          spat_bias = TRUE,
          samp_weights = NULL,
          bias_level = 3,
          rep_visit = 2.5,
          maxVisit = 10,
          reli_score = as.character(exp_cond[i, 3])
        )
      
      # 'Real' species richness
      real_sr <- spp_occ %>%
        dplyr::left_join(hexRelate, . , by = 'id') %>%
        dplyr::select(c('id2', grep(names(.), pattern = 'pa_'))) %>%
        dplyr::group_by(id2) %>%
        dplyr::summarise(across(starts_with("pa_"), ~ max(.))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(sr = sum(across(starts_with("pa_")))) %>%
        dplyr::select(c('id2', 'sr'))
      
      # Loop through cells calculating coverage
      surv_ex <- surv_ex %>%
        dplyr::left_join(hexRelate, . , by = 'id') %>%
        dplyr::mutate(id2 = as.numeric(id2))  %>%
        dplyr::as_tibble()
      
      # Results
      res_out <- data.frame(
        id2 =  unique(surv_ex$id2),
        n_spp = exp_cond[i, 1],
        p_sites = exp_cond[i, 2],
        reli_score = as.character(exp_cond[i, 3]),
        n_samp = NA,
        real_sr = NA,
        obs_spp = NA,
        chao = NA,
        chao.se = NA,
        jack1 = NA,
        jack1.se = NA,
        jack2 = NA,
        boot = NA,
        boot.se = NA
      )
      
      # Loop through coarser cells
      for (x in 1:nrow(res_out)) {
        # subset cell
        tmp <- surv_ex %>%
          dplyr::as_tibble() %>%
          dplyr::filter(id2 == res_out[x, 1]) %>%
          dplyr::select(grep(names(.), pattern = 'det_'))  %>%
          replace(is.na(.), 0)
        
        # Rarefaction
        #rareHill <- rarefaction.sample(tmp, q = 0, method = "coverage")
        specRes <- specpool(tmp, smallsample = T)
        
        # Output
        res_out[x, 5] <- specRes$n
        res_out[x, 6] <- real_sr %>%
          dplyr::filter(id2 == res_out[x, 1]) %>%
          dplyr::select('sr')
        res_out[x, 7] <- specRes$Species
        res_out[x, 8] <- specRes$chao
        res_out[x, 9] <- specRes$chao.se
        res_out[x, 10] <- specRes$jack1
        res_out[x, 11] <- specRes$jack1.se
        res_out[x, 12] <- specRes$jack2
        res_out[x, 13] <- specRes$boot
        res_out[x, 14] <- specRes$boot.se
        
      }
      
      return(res_out)
      
    })
    
    res_out <- do.call(rbind, res_out)
    
  }

#' Calculate first derivatives of records x time
#' 
#' @param nRec A dataframe containing the number of records per year
#' @param spp Names of the speices group.
#' 
#' @return A dataframe showing the periods of most rapid increase in records.
firstDeriv_nRec <- function(nRec, spp) { 
  
  DF <- nRec %>%
    filter(sppGrp == spp) %>%
    mutate(Yr = as.numeric(Yr))
  
  # Fit poisson model
  mod <-
    gam(log10(nRec) ~ s(Yr), data = DF, family = 'gaussian')
  plot(mod)
  
  # new data for prediction
  newdf <- with(DF, data.frame(Yr = unique(Yr)))
  
  # Predict to new data and convert to the response scale
  ilink <- family(mod)$linkinv      # inverse link function
  DF2 <- predict(mod, newdata = newdf, se.fit = TRUE)
  DF2 <- with(DF2,
              data.frame(
                newdf,
                response = ilink(fit),
                lwr = ilink(fit - 2 * se.fit),
                upr = ilink(fit + 2 * se.fit)
              ))
  
  # Calculate first derivatives
  dv <- derivatives(mod, newdata = newdf, type = "central")
  dv <- dv[, c('derivative', 'lower', 'upper')]
  DF2 <- cbind(DF2, dv)
  
  # Identify significant positive first derivatives
  DF2$fdYrs <- 0
  DF2$fdYrs[which(DF2$derivative > 0 & DF2$lower > 0)] <- 1
  
  # Identify clusters of growth
  gr <-
    data.frame(
      start = 1,
      end = cumsum(rle(DF2$fdYrs)$length),
      v = rle(DF2$fdYrs)$values
    )
  for (i in 2:nrow(gr))
    gr[(i), 1] <- gr[(i - 1), 2] + 1
  gr <- gr[gr$v == 1,]
  gr$id <- cumsum(gr$v)
  for (i in 1:nrow(gr))
    DF2$fdYrs[gr$start[i]:gr$end[i]] <- gr$id[i]
  DF2$sppGrp <- spp
  DF2
  
}

#' Run cover estimates for each species group and grain
#' 
#' @param sppLookUp Dataframe with species group names and names for plotting.
#' @param yearStart Numeric giving start year for analysis.
#' @param grains A numeric or numeric vector giving grains (area) to run.
#' @param yearType Either continuous or decadal.
#' @param outPath A file path for writing out the results.
#' @param collateOnly If the models have been run but the results need to be 
#' collated this function can skip the model running phase by setting to TRUE.
#' 
#' @return NULL
run_cover <- function(sppLookUp = sppNameTab,
                      yearStart = 1960,
                      grains = c(1, 4, 16, 64),
                      yearType = c("continuous", "decadal"),
                      outPath = 'Outputs/Temporal_coverage_R1/',
                      collateOnly = FALSE) {
  sppGrps <- unique(sppLookUp[, 2])
  
  if (!collateOnly) {
    for (sppGrpN in sppGrps) {
      #~# Subset to single species group
      sppSelGrp <- NULL
      erccis_sppGrp <-
        species.group.select(erccis_data,
                             sppLookUp[, 1][which(sppLookUp[, 2] == sppGrpN)],
                             cornwall_sf,
                             "Cornwall",
                             sppSelect = sppSelGrp)
      #~# Run for each grain
      for (grn in grains) {
        hexGrain <- makeHexGrid(cornwall_sf, grn)
        covEstOut <- estCovSppGrp2(
          erccis_sppGrp,
          sppGrpN,
          timeSlice = 1,
          yearStart = yearStart,
          region_sf = cornwall_sf,
          grain = grn,
          hexGrid = hexGrain,
          yearType = yearType,
          outPath = outPath
        )
      }
    }
  }
  for (i in 1:length(sppGrps)) {
    sppName <- as.character(sppGrps[i])
    sppList <-
      list.files(outPath, pattern = paste0('^Coverage_', sppName))
    for (j in c(1, 4, 16, 64)) {
      grnSppList <-
        sppList[grep(sppList, pattern = paste0('grain_', j, '_'))]
      sppCovAll <- lapply(grnSppList, function(x) {
        tmp_file <- read.csv(paste0(outPath, x))
        tmp_file$sppGrp <- sppName
        tmp_file$grain <- j
        tmp_file
      })
      sppCovAll <- do.call(bind_rows, sppCovAll)
      op <-
        paste0(outPath, 'Temporal_coverage_', sppName, '_', j, '.csv')
      write.csv(sppCovAll,  file = op)
    }
  }
}

#' Estimating sampling coverage using the iNext functions
#' 
#' Note: Function runs in parallel using snowfall
#' 
#' @param erccis_sppGrp ERCCIS data filtered for a particular species group.
#' @param sppGrp Name of the species group.
#' @param timeSlice A numeric value specifying the length on each timeslice 
#' (in years).
#' @param yearStart Numeric giving start year for analysis.
#' @param region_sf An sf object defining the spatial region.
#' @param grain Numeric giving the grain for the analysis.
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param yearType Either continuous or decadal.
#' @param outPath A file path for writing out the results.
#' 
#' @return NULL
estCovSppGrp2 <- function(erccis_sppGrp,
                          sppGrp,
                          timeSlice = 1,
                          yearStart = 1970,
                          region_sf,
                          grain = 1,
                          hexGrid = NULL,
                          yearType = c("continuous", "decadal"),
                          outPath) {
  
  #~# Select all years
  siteXSppYr <- siteXSpeciesXTime(erccis_sppGrp, 1, 
                                  region_sf, 
                                  grain, 
                                  NULL, 
                                  hexGrid)
  siteXSppYr$site <- sapply(siteXSppYr$siteTime, 
                            function(x) strsplit(x, split = "_")[[1]][1])
  siteXSppYr$year <- sapply(siteXSppYr$siteTime, 
                            function(x) strsplit(x, split = "_")[[1]][2])

  #- Set up snowfall
  sfInit( parallel = TRUE, cpus = 7 )
  sfLibrary(gtools)
  sfLibrary(tidyverse)
  sfLibrary(iNEXT.4steps)
  sfExport(list = c("siteXSppYr","timeSlice","yearStart",
                    "outPath","sppGrp","grain","yearType"))
  
  #~# Species coverage estimation yr <- 1971
  res_all <- sfLapply(seq((yearStart + timeSlice),
                          2019, timeSlice), function(yr) {
                            print(yr)
                            
                            fname <- paste0(outPath,
                                            "Coverage_",
                                            sppGrp,
                                            "_ts_",
                                            timeSlice,
                                            "_grain_",
                                            grain,
                                            "_yr_",
                                            yr,
                                            ".csv")
                            
    if(!file.exists(fname)) {
      
      #~# Sample increasing amounts of time
      if ( yearType == "continuous" ) {
        ts_dat <- siteXSppYr[siteXSppYr$year <= yr, ]
      }
      if ( yearType == "decadal" ) {
        ts_dat <- siteXSppYr[siteXSppYr$year %in% seq(yr - 10, yr), ]
      }
      
      #~# Calculate site based sampling coverage i = 13
      covRes <- data.frame(
        site = mixedsort(unique(siteXSppYr$site)),
        year = yr,
        n_samp = NA,
        obs_spp = NA,
        SC_q0 = NA,
        SC_q1 = NA,
        SC_q2 = NA,
        sr_asym_est = NA,
        sr_asym_se = NA,
        sr_asym_LCL = NA,
        sr_asym_UCL = NA,
        sr_nonasym_q0 = NA,
        sr_nonasym_q1 = NA,
        sr_nonasym_q2 = NA
      )
      covRes$site <- as.numeric(as.character(covRes$site))
      
      #~# Run through all sites
      for(i in 1:nrow(covRes)) {
        
        # Subset to time period
        ts_dat_tmp <- ts_dat %>%
          filter(site == covRes[i,1]) %>%
          select(!c('siteTime','site','year'))
        ts_dat_tmp[ts_dat_tmp > 1] <- 1
        
        # N species
        covRes[i,"obs_spp"] <- length(which(colSums(ts_dat_tmp) > 0))
        covRes[i,"n_samp"] <- nrow(ts_dat_tmp)
        
        if( covRes[i,"n_samp"] >= 10 & covRes[i,"obs_spp"] > 1 ) {
          
          #specRes <- specpool(ts_dat_tmp, smallsample = TRUE)
          ts_freq_lst <- c(nrow(ts_dat_tmp),as.numeric(colSums(ts_dat_tmp)))
          
          possibleError1 <- 
            tryCatch(qProf <- iNEXT4steps(ts_freq_lst, 
                                          diversity = "TD",
                                          datatype = "incidence_freq"),
                     error=function(e) e) 
          
          if(!inherits(possibleError1, "error")){
            
            # Sample completeness profiles
            covRes[i,"SC_q0"] <- 
              qProf$summary$`STEP1. Sample completeness profiles`$`q = 0`
            covRes[i,"SC_q1"] <- 
              qProf$summary$`STEP1. Sample completeness profiles`$`q = 1`
            covRes[i,"SC_q2"] <- 
              qProf$summary$`STEP1. Sample completeness profiles`$`q = 2`
            
            # Asymptotic analysis species richness
            covRes[i,"sr_asym_est"] <- 
              qProf$summary$`STEP2. Asymptotic analysis`$Estimator[1]
            covRes[i,"sr_asym_se"] <- 
              qProf$summary$`STEP2. Asymptotic analysis`$s.e.[1]
            covRes[i,"sr_asym_LCL"] <- 
              qProf$summary$`STEP2. Asymptotic analysis`$LCL[1]
            covRes[i,"sr_asym_UCL"] <- 
              qProf$summary$`STEP2. Asymptotic analysis`$UCL[1]
            
            # Non-asymptotic coverage-based rarefaction & extrapolation analysis
            covRes[i,"sr_nonasym_q0"] <- qProf$summary[[3]]$`q = 0`
            covRes[i,"sr_nonasym_q1"] <- qProf$summary[[3]]$`q = 1`
            covRes[i,"sr_nonasym_q2"] <- qProf$summary[[3]]$`q = 2`
            
          }
          
        }
        
      }
      
      # Save each year (otherwise it bogs down)
      write.csv(covRes, fname, row.names = FALSE)
      
    }
    
    return(NULL)
    
  })
  
  sfStop()
  
  return(NULL)
  
}

#' C-index of coverage
#' 
#' @param filePath The file path to the folder containing the coverage estimates
#' @param q The order of q, which determines the sensitivity of coverage 
#' estimates to infrequent species, with q = 1 giving the sample coverage (i.e.
#' the proportion of the total number of individuals in the entire assemblage 
#' that belong to detecte species).
#' 
#' @return
cIndexEstCov2 <- function(filePath,
                          q = "SC_q0") {
  #~# Extract all species data
  f_all <-
    list.files(filePath, pattern = '^Temporal_coverage', full.name = TRUE)
  sampCov_df <- lapply(f_all, function(f) {
    print(f)
    dat <- read.csv(f)[, -1]
  })
  sampCov_df <- do.call(rbind, sampCov_df)
  
  #~# total cells
  nCell <- sampCov_df %>%
    select(site, grain) %>%
    group_by(grain) %>%
    distinct() %>%
    group_by(grain) %>%
    tally()
  colnames(nCell)[2] <- 'N'
  
  #~# add total cells
  sampCov_df <- left_join(sampCov_df, nCell, by = 'grain')
  
  #~# Apply threshold sample sizes
  sampCov_df[[q]][sampCov_df$grain == 1 &
                    sampCov_df$n_samp  < 10] <- NA
  sampCov_df[[q]][sampCov_df$grain == 4 &
                    sampCov_df$n_samp  < 25] <- NA
  sampCov_df[[q]][sampCov_df$grain == 16 &
                    sampCov_df$n_samp  < 50] <- NA
  sampCov_df[[q]][sampCov_df$grain == 64 &
                    sampCov_df$n_samp  < 50] <- NA
  
  #~# Select q
  sampCov_df <-
    sampCov_df[, c("site", "year", "sppGrp", q, "grain", "N", "n_samp")]
  names(sampCov_df)[4] <- "estCov"
  
  #sampCov_df %>% filter(sppGrp == "Birds", year == "2019", grain == 64)
  
  #~# H-index i <- 0.9
  cIndex <- lapply(seq(0, 1, 0.01), function(i) {
    covRes <- sampCov_df %>%
      mutate(estCov = round(estCov, 2)) %>%
      group_by(year, sppGrp, grain) %>%
      filter(estCov > i) %>%
      group_by(year, sppGrp, grain, N) %>%
      summarise(cInd = n(), .groups = 'drop') %>%
      mutate(cInd = cInd / N) %>%
      group_by(sppGrp, grain) %>%
      filter(cInd > i) %>%
      group_by(sppGrp, grain) %>%
      filter(cInd == min(cInd))
  })
  cIndex <- do.call(rbind, cIndex)
  cIndex$grain <- as.factor(cIndex$grain)
  cIndex <- cIndex %>%
    distinct() %>%
    arrange(sppGrp, grain, year) %>%
    group_by(sppGrp, grain, year) %>%
    filter(cInd == max(cInd))
  cIndex$sppGrp <- sub('_', ' ', cIndex$sppGrp)
  cIndex$sppGrp <- factor(cIndex$sppGrp,
                          levels = sub("_", " ", unique(sppNameTab[, 2])))
  
  cIndex
}

#' Completeness index panel plots
#'
#' @param filePath The file path to the folder containing the coverage 
#' estimates.
#' @param outPath The destination and file name for plots and tables.
#' @param q The order of q, which determines the sensitivity of coverage 
#' estimates to infrequent species, with q = 1 giving the sample coverage (i.e.
#' the proportion of the total number of individuals in the entire assemblage 
#' that belong to detecte species).
#'
ci_panel_plot <- function(filePath,
                          outPath,
                          q = c('SC_q0','SC_q1','SC_q2')) {
                            
  ci_dat <- cIndexEstCov2(filePath, q)
  
  #- Create panel plot
  fd_ci_dat <- run_firstDerivCoverage(ci_dat)
  fd_ci_dat$sppGrp <- factor(fd_ci_dat$sppGrp,
                             levels = sub("_", " ", unique(sppNameTab[, 2])))
  
  
  ci_panel <- ggplot(ci_dat) +
    geom_line(
      aes(
        x = year,
        y = response,
        group = interaction(fdYrs, grain)
      ),
      colour = 'Red',
      size = 2.5,
      alpha = 0.4,
      data = fd_ci_dat[which(fd_ci_dat$fdYrs > 0), ]
    ) +
    geom_point(aes(
      x = year,
      y = cInd,
      group = grain,
      colour = grain
    ),
    size = 1,
    alpha = 0.3) +
    geom_line(
      aes(
        x = year,
        y = response,
        group = grain,
        colour = grain
      ),
      size = 0.5,
      data = fd_ci_dat
    ) +
    scale_color_colorblind(name = expression('Grain (km' ^ 2 * ')')) +
    scale_linetype(name = expression('Grain (km' ^ 2 * ')')) +
    theme_bw(base_size = 16) %+replace%
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45)) +
    ylim(0, 1) + xlim(1960, 2020) + ylab('Completeness Index') + xlab('Year') +
    guides(
      colour = guide_legend(override.aes = list(size = 1)),
      linetype = guide_legend(override.aes = list(size = 1))
    ) +
    facet_wrap( ~ sppGrp, ncol = 5)
  
  ggsave(
    filename = paste0(outPath, ".png"),
    plot = ci_panel,
    width = 12,
    height = 8
  )
  
  #- Tabulated final outcomes
  ci_table <- ci_dat %>%
    drop_na() %>%
    group_by(sppGrp, grain) %>%
    filter(year == max(year)) %>%
    select(sppGrp, grain, cInd) %>%
    arrange(sppGrp, grain) %>%
    pivot_wider(id_cols = 'sppGrp',
                names_from = 'grain',
                values_from = 'cInd')
  write.csv(ci_table,
            file = paste0(outPath, ".csv"),
            row.names = FALSE)
  
}

#' Proportion of cells exceeding threshold for sampling completeness
#' 
#' @param filePath The file path to the folder containing the coverage 
#' estimates.
#' @param threshold The threshold for sampling completeness.
#' @param q The order of q, which determines the sensitivity of coverage
#' estimates to infrequent species, with q = 1 giving the sample coverage (i.e.
#' the proportion of the total number of individuals in the entire assemblage
#' that belong to detected species).
#'
#' @import
#'
#' @return 
propCovThres2 <- function(filePath = 'Outputs/Temporal_coverage_R1_continuous/', 
                          threshold = 0.8,
                          q = 'SC_q0',
                          max_pc = FALSE) {
  
  #~# Extract all species data
  f_all <-
    list.files(filePath, pattern = '^Temporal_coverage', full.name = TRUE)
  sampCov_df <- lapply(f_all, function(f) {
    dat <- read.csv(f)[, -1]
  })
  sampCov_df <- do.call(rbind, sampCov_df)
  
  #~# total cells
  nCell <- sampCov_df %>%
    select(site, grain) %>%
    distinct() %>%
    group_by(grain) %>%
    summarise(N = n())
  
  #~# add total cells
  sampCov_df <- sampCov_df %>%
    left_join(nCell, by = 'grain')
  
  #~# Apply threshold sample sizes
  sampCov_df[[q]][sampCov_df$grain == 1 &
                    sampCov_df$n_samp  < 10] <- NA
  sampCov_df[[q]][sampCov_df$grain == 4 &
                    sampCov_df$n_samp  < 25] <- NA
  sampCov_df[[q]][sampCov_df$grain == 16 &
                    sampCov_df$n_samp  < 50] <- NA
  sampCov_df[[q]][sampCov_df$grain == 64 &
                    sampCov_df$n_samp  < 50] <- NA
  
  #~# For expanding all years
  expYears <- expand.grid(
    year = seq(1960, 2019, 1),
    sppGrp = unique(sppNameTab$SpeciesGroupFile),
    grain = c(1, 4, 16, 64)
  )
  
  #~# Select q
  sampCov_df <-
    sampCov_df[, c("site", "year", "sppGrp", q, "grain", "N", "n_samp")]
  names(sampCov_df)[4] <- "estCov"
  
  if (!max_pc) {
    #~# Calculate year threshold coverage is exceeded
    covThresYr <- sampCov_df %>%
      #~# Select columns
      select(sppGrp, grain, year, site, estCov, N) %>%
      #~# Remove years/cells with no cover estimated
      filter(!is.na(estCov)) %>%
      arrange(sppGrp, grain, site, year) %>%
      #~# Filter cells/years that cross threshold
      filter(estCov >= threshold) %>%
      #~# Number of instances per year
      group_by(sppGrp, grain, year) %>%
      summarise(n = n(),
                pc = n / N,
                .groups = 'drop') %>%
      distinct() %>%
      #~# Expand to include all years
      full_join(expYears, covThresYr, by = c('year', 'sppGrp', 'grain')) %>%
      #~# Arrange by year
      arrange(sppGrp, grain, year)
  
    #~# Fill in NAs with preceding value
    covThresYr$pc[covThresYr$year == 1960] <- 0
    for (i in 2:nrow(covThresYr))
      if (is.na(covThresYr$pc[i]))
        covThresYr$pc[i] <- covThresYr$pc[i - 1]
    covThresYr$grain <- as.factor(covThresYr$grain)
    covThresYr$sppGrp <- sub('_', ' ', covThresYr$sppGrp)
    covThresYr$sppGrp <- factor(covThresYr$sppGrp,
                                levels = sub("_", " ", unique(sppNameTab[, 2])))
  } else {
    covThresYr <- sampCov_df
  }
    
  #~# Return for plotting
  covThresYr
  
}

#' Compiles results and produces panel and table for the proportion of cells
#'  approaching sampling completeness (e.g. >80%)
#' 
#' @param filePath The file path to the folder containing the coverage 
#' estimates.
#' @param outPath The destination and file name for plots and tables.
#' @param threshold The threshold for sampling completeness.
#' @param q The order of q, which determines the sensitivity of coverage
#' estimates to infrequent species, with q = 1 giving the sample coverage (i.e.
#' the proportion of the total number of individuals in the entire assemblage
#' that belong to detected species).
#'
#' @export NULL
pc_panel_plot <- function(filePath,
                          outPath,
                          threshold = 0.8,
                          q = 'SC_q0') {
  pc_dat <-
    propCovThres2(filePath,
                  threshold,
                  q)
  fd_pc_dat <- run_firstDerivCoverage(pc_dat)
  fd_pc_dat$sppGrp <- factor(fd_pc_dat$sppGrp,
                              levels = sub("_", " ", unique(sppNameTab[, 2])))
  
  pc_cont_p <- ggplot(pc_dat) +
    geom_line(
      aes(
        x = year,
        y = response,
        group = interaction(fdYrs, grain)
      ),
      colour = 'Red',
      size = 2.5,
      alpha = 0.4,
      data = fd_pc_dat[which(fd_pc_dat$fdYrs > 0), ]
    ) +
    geom_point(
      aes(
        x = year,
        y = pc,
        group = as.factor(grain),
        colour = as.factor(grain)
      ),
      size = .75,
      alpha = 0.5
    ) +
    geom_line(
      aes(
        x = year,
        y = response,
        group = grain,
        colour = grain
      ),
      size = 0.5,
      data = fd_pc_dat
    ) +
    scale_color_colorblind(name = expression('Grain (km' ^ 2 * ')')) +
    scale_linetype(name = expression('Grain (km' ^ 2 * ')')) +
    theme_bw(base_size = 16) %+replace%
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45)) +
    ylim(0, 1) + xlim(1960, 2020) +
    ylab('Proportion cells >80% completeness') + xlab('Year') +
    guides(colour = guide_legend(override.aes = list(size = 1)),
           linetype = guide_legend(override.aes = list(size = 1))) +
    facet_wrap( ~ sppGrp, ncol = 5, scales = "fixed")
 
  ggsave(
    filename = paste0(outPath, ".png"),
    plot = pc_cont_p,
    width = 12,
    height = 8
  )
  
  #-- Create table
  pc_table <- pc_dat %>%
    drop_na() %>%
    group_by(sppGrp, grain) %>%
    filter(pc == max(pc)) %>%
    select(sppGrp, grain, pc) %>%
    distinct() %>%
    arrange(sppGrp, grain) %>%
    pivot_wider(id_cols = 'sppGrp',
                names_from = 'grain',
                values_from = 'pc')
  write.csv(pc_table,
            file = paste0(outPath, ".csv"),
            row.names = FALSE)
  
  return(NULL) 
  
}

#' Calculate first derivatives of coverage 
#' 
#' @param covDat A dataframe with coverage estimates x time.
#' @param spp The name of the species group
#' @param grn The grain (area) used in modelling.
#' 
#' @return A dataframe showing the periods of most rapid increase in coverage.
firstDerivCoverage <- function(covDat, spp, grn) {
  print(spp)
  # Prepare data
  yearsAll <- data.frame(year = 1961:2019)
  DF <- covDat %>%
    filter(sppGrp == spp & grain == grn)
  DF <- left_join(yearsAll, DF, by = 'year')
  names(DF)[which(names(DF) %in% c('cInd', 'pc'))] <- 'value'
  if (is.na(DF$value[1]))
    DF$value[1] <- 0
  for (i in 2:nrow(DF))
    if (is.na(DF$value[i]))
      DF$value[i] <- DF$value[(i - 1)]
  DF <- DF %>% select(year, sppGrp, grain, value)
  
  if (sum(DF$value) > 0) {
    # Model using betar gam
    mod <-
      gam(value ~ s(year), data = DF, family = betar())
    plot(mod)
    # New data for prediction
    newdf <- with(DF, data.frame(year = unique(year)))
    # Predict to new data and convert to the response scale
    ilink <- family(mod)$linkinv      # inverse link function
    DF2 <- predict(mod, newdata = newdf, se.fit = TRUE)
    DF2 <- with(DF2,
                data.frame(
                  newdf,
                  response = ilink(fit),
                  lwr = ilink(fit - 2 * se.fit),
                  upr = ilink(fit + 2 * se.fit)
                ))
    # Calculate first derivatives
    DF2 <- cbind(DF2,
                 derivatives(mod,
                             newdata = newdf,
                             type = "central")[, c('derivative', 'lower', 'upper')])
    # Identify significant positive first derivatives
    DF2$fdYrs <- 0
    DF2$fdYrs[which(DF2$derivative > 0 & DF2$lower > 0)] <- 1
    DF2$fdYrs[which(DF2$derivative <
                      quantile(DF2$derivative[which(DF2$fdYrs == 1)],
                               prob = 0.5) &  DF2$fdYrs == 1)] <- 0
    
    # Identify clusters of growth
    gr <- data.frame(
      start = 1,
      end = cumsum(rle(DF2$fdYrs)$length),
      v = rle(DF2$fdYrs)$values
    )
    if (nrow(gr) > 1) {
      for (i in 2:nrow(gr))
        gr[(i), 1] <- gr[(i - 1), 2] + 1
      gr <- gr[gr$v == 1, ]
      gr$id <- cumsum(gr$v)
      for (i in 1:nrow(gr))
        DF2$fdYrs[gr$start[i]:gr$end[i]] <- gr$id[i]
    }
    DF2$sppGrp <- spp
    DF2$grain <- as.factor(grn)
    DF2
  }
  
}

#' Run the derivative function.
#' 
#' @param sppNames Vector of species names that link to covDat names
#' 
#' @return A dataframe showing the periods of most rapid increase in coverage.
run_firstDerivCoverage <-
  function(covDat, sppNames = unique(sppNameTab[, 2])) {
    fdSpAll <- lapply(sppNames, function(spp) {
      spp <- sub('_', ' ', spp)
      spGrn <- lapply(c(1, 4, 16, 64), function(grn) {
        fd <- firstDerivCoverage(covDat, spp, grn)
      })
      spGrn <- do.call(bind_rows, spGrn)
      spGrn
    })
    fdSpAll <- as_tibble(do.call(bind_rows, fdSpAll))
  }

#' Proportion of cells exceeding threshold per cell and species group
#' 
#' @param sppLookUp Dataframe with species group names and names for plotting.
#' @param grain The grain (area) used in modelling.
#' @param threshold The threshold for sampling completeness.
#' 
#' @return sf object showing the inventory completeness per cell x species grp.
propCovThresConsensus <- function(sppNameLookUp,
                                  grain = 4,
                                  threshold = 0.8) {
  
  #~# File path combinations
  covFiles <- expand.grid(sppGrp = unique(sppNameLookUp$SpeciesGroupFile), 
                          grain = grain)
  
  #~# Extract all species data
  sampCov_df <- lapply(1:nrow(covFiles), function(i) {
    dat <- read.csv(paste0('Outputs/Coverage_',
                           covFiles[i,1],'_Grain_',
                           covFiles[i,2],'.csv'))[,-1]
  })
  sampCov_df <- do.call(rbind,sampCov_df) 
  
  #~# total cells
  nCell <- sampCov_df %>%
    select(site, grain) %>%
    group_by(grain) %>%
    filter(site == max(site)) %>% 
    distinct()
  colnames(nCell)[1] <- 'N'
  
  #~# add total cells
  sampCov_df <- sampCov_df %>%
    left_join(nCell, by = 'grain')
  
  #~# For expanding all years
  expYears <- expand.grid(year = seq(1960,2020,1),
                          sppGrp = unique(sppNameLookUp$SpeciesGroupFile),
                          grain = grain)
  
  #~# Calculate year threshold coverage is exceeded
  covThresYr <- sampCov_df %>%
    #~# Select columns
    select(sppGrp, year, site, estCov, N) %>%
    #~# Remove years/cells with no cover estimated
    filter(!is.na(estCov)) %>%
    # #~# Sometimes coverage appears high in the first few years because of small samples
    group_by(sppGrp, site) %>%
    filter(estCov <= last(estCov)) %>% 
    #~# Filter cells/years that cross threshold
    mutate(estCov = ifelse(estCov > threshold, 1, 0)) %>%
    #~# ... and then get the first year this occurred
    filter(year == max(year))  %>% 
    #~# Arrange by year
    arrange(sppGrp, site) %>% 
    #~# Remove unnecessary columns
    select(sppGrp, site, estCov) %>%
    #~# Get only cells with coverage
    filter(estCov == 1) %>%
    #~# Summarise
    group_by(site) %>%
    summarise(nCov = sum(estCov)) %>%
    rename(id = site) %>%
    mutate(id = as.character(id))
  
  #~# Hex grid
  hexGrn <- st_make_grid(cornwallbound_sf, 
                         cellsize = (1074.57 * sqrt(grain)), 
                         square = FALSE) %>%
    st_as_sf() %>%  mutate(id = row.names(.), 
                           grain = grain) # plot(hexGrid_256) get(paste0('hexGrid_',grn))
  
  #~# Spatial join
  covThresYrHex <- left_join(hexGrn, covThresYr, by = 'id')
  covThresYrHex$nCov[is.na(covThresYrHex$nCov)] <- 0
  covThresYrHex
  
}

#' Maximum coverage estimate
#' 
#' @param filePath File path of coverage estimates.
#' @param yearMin The minimum year to check across to find the max. coverage.
#' @param q The order of q, which determines the sensitivity of coverage
#' estimates to infrequent species, with q = 1 giving the sample coverage (i.e.
#' the proportion of the total number of individuals in the entire assemblage
#' that belong to detected species).
#' 
#' @return A dataframe with max coverage estimate per cell.
maxEstCov <-
  function(filePath = 'Outputs/Temporal_coverage_R1_continuous/',
           yearMin = 1960,
           q = "SC_q1") {
    #~# Extract all species data
    f_all <-
      list.files(filePath, pattern = '^Temporal_coverage', full.name = TRUE)
    sampCov_df <- lapply(f_all, function(f) {
      dat <- read.csv(f)[, -1]
    })
    sampCov_df <- do.call(rbind, sampCov_df)
    
    #~# Apply threshold sample sizes
    sampCov_df[[q]][sampCov_df$grain == 1 &
                      sampCov_df$n_samp  < 10] <- NA
    sampCov_df[[q]][sampCov_df$grain == 4 &
                      sampCov_df$n_samp  < 25] <- NA
    sampCov_df[[q]][sampCov_df$grain == 16 &
                      sampCov_df$n_samp  < 50] <- NA
    sampCov_df[[q]][sampCov_df$grain == 64 &
                      sampCov_df$n_samp  < 50] <- NA
    
    #~# Select q
    sampCov_df <-
      sampCov_df[, c("site", "year", "sppGrp", q, "grain", "n_samp", "obs_spp")]
    names(sampCov_df)[4] <- "estCov"
    
    sampCov_df$estCov[is.na(sampCov_df$estCov)] <- 0
    sampCov_df <- sampCov_df  %>%
      filter(year >= yearMin)  %>%
      select(site, grain, obs_spp, estCov, n_samp, year, sppGrp) %>%
      group_by(site, sppGrp, grain) %>%
      summarise(covMax = max(estCov, na.rm = TRUE),
                nRec = max(obs_spp, na.rm = TRUE))
    
    #~# All combination of site and grains
    allCellGrn <- sampCov_df %>%
      ungroup() %>%
      select(site, grain) %>%
      distinct()
    sppGrn <-
      expand.grid(grain = c(1, 4, 16, 64),
                  sppGrp = unique(sppNameTab[, 2]))
    allCellGrn <- left_join(allCellGrn, sppGrn, by = 'grain')
    #allCellGrn %>% group_by(grain, sppGrp) %>% tally()
    
    #~# Combine and fill out
    sampCov_df <-
      left_join(allCellGrn, sampCov_df, by = c('site', 'grain', 'sppGrp'))
    #sampCov_df %>% group_by(sppGrp, grain) %>% tally()
    
    #~# Return
    sampCov_df
    
  }

#~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~#
#~# #~# #~# ENVIRONMENTAL DATA PROCESSING FUNCTIONS #~# #~# #~# 

#' Distance decay to environmental features
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import sf
#'
#' @return Dataframe
distance.decay <- function(enviro_sf,
                           region_sf,
                           hexCellArea,
                           hexGrid = NULL,
                           var_name)
{
  #-- Create hexagonal grid with id attribute
  if (is.null(hexGrid))
    hexGrid <- makeHexGrid(region_sf, hexCellArea)
  #-- Clip to region
  enviro_clip <- st_crop(enviro_sf, region_sf)
  #-- Calculate distance from each grid cell to poly / point
  enviro_dist <- st_distance(hexGrid, enviro_clip)
  enviro_hex_dist <-
    cbind(hexGrid, value = apply(enviro_dist, 1, min) / 1000)
  enviro_hex_dist$grain <- hexCellArea
  #-- Drop geometry and return
  enviro_hex_dist <- st_drop_geometry(enviro_hex_dist) %>%
    arrange(id, grain, value) %>%
    select(id, value, grain)
  enviro_hex_dist$var_name <- var_name
  enviro_hex_dist
  
}

#' Count feature occurrence within cells
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import sf
#'
#' @return An sf object.
#'
#' @export Dataframe
feature.occurence <- function(enviro_sf,
                              region_sf,
                              hexCellArea,
                              hexGrid = NULL,
                              var_name)
{
  #~# Create hexagonal grid with id attribute
  if (is.null(hexGrid))
    hexGrid <- makeHexGrid(region_sf, hexCellArea)
  #~# Count features within cells
  featureCount <-
    st_join(hexGrid, enviro_sf, join = st_intersects) %>%
    group_by(id, grain) %>%
    summarise(value = as.numeric(n()))
  featureCount$grain <- hexCellArea
  #~# Print test plot
  # print(ggplot(featureCount)+geom_sf(aes(fill=value))+
  #         scale_fill_distiller(palette = 'Greys', direction = 1))
  #~# Prep and return
  featureCount <- featureCount %>%
    st_drop_geometry() %>%
    arrange(id, grain, value)
  featureCount$var_name <- var_name
  featureCount
  
}

#~# Count feature occurrence within cells enviro_sf <- highway; region_sf <- cornwallbound_sf; hexCellArea <- 4; hexGrid <- hex4
#' 
#' @param 
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' 
#' @import 
#' 
#' @return 
feature.length <- function(enviro_sf,
                           region_sf,
                           hexCellArea,
                           hexGrid = NULL,
                           var_name)
{
  #~# Create hexagonal grid with id attribute
  if (is.null(hexGrid))
    hexGrid <- makeHexGrid(region_sf, hexCellArea)
  
  #~# Count features within cells
  fl <- st_join(enviro_sf, hexGrid, join = st_intersects)
  
  #~# Calculate road/path lengths
  fl.res <- data.frame(id = na.omit(unique(hexGrid$id)),  value = NA)
  for (i in fl.res$id) {
    flg <- fl %>% filter(id == i)
    lKm <-
      as.numeric(sum(st_length(st_cast(
        flg, 'LINESTRING'
      ))) / 1000)
    fl.res[which(fl.res[, 1] == i), 2] <- lKm
  }
  fl.res$grain <- hexCellArea
  fl.res <- fl.res %>% arrange(id, grain, value)
  fl.res$var_name <- var_name
  fl.res
  
}

#' Count feature occurrence within cells
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import  sf
#'
#' @return An sf object.
#'
#' @export Dataframe
feature.cover <- function(enviro_sf, 
                          region_sf, 
                          hexCellArea,
                          hexGrid = NULL,
                          var_name) 
{
  
  #--Create hexagonal grid with id attribute
  if(is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf,hexCellArea)
  
  #-- Count features within cells
  fa <- st_intersection(hexGrid, enviro_sf) 
  fa.out <- data.frame(id = unique(hexGrid$id), grain = hexCellArea, value = NA)

  #-- run through cells
  for(i in unique(hexGrid$id)){
    print(i)
    fai <- fa %>%
      filter(id == i) 
    if(nrow(fai)!=0){
      fai <- fai %>%
        st_cast() %>%
        st_area()
      fa.out[which(fa.out[,1]==i),3] <- round(as.numeric(fai / 10^6),3)
    }
  }
  fa.out$var_name <- var_name
  fa.out
  
}

#~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~# #~#
#~# #~# #~#            Spatial modelling            #~# #~# #~#

#' Prepare data for modelling
prepareModelDat <- function(hexGrid, grn, covEnvDat) {
  
  # Get coordinates
  XYid <- data.frame(st_coordinates(st_centroid(hexGrid)), 
                     id = as.numeric(hexGrid$id))
  
  #~# Process data to grain and threshold
  hexGrainArea <- hexGrid %>%
    st_drop_geometry() %>%
    mutate(id = as.numeric(id), area = as.numeric(area)) %>%
    dplyr::select('id', 'grain', 'X', 'Y', 'area')
  
  # Convert to 0 and 1 at a threshold
  maxEnvGrn <- covEnvDat %>%
    group_by(sppGrp, grain, id) %>%
    mutate(covComp = as.integer(ifelse(covMax > 0.8, 1, 0))) %>%
    filter(!is.na(sppGrp), grain == grn) %>% 
    full_join(hexGrainArea, by = c('id', 'grain')) %>%
    replace_na(list(covComp = 0, highway = 0, nSchemes = 0)) %>%
    filter(area > 0) %>%
    mutate(area = area / {{grn}})
  
  # Return
  maxEnvGrn
  
}

#' Create neighbourhood weighting object
create.nb.mat <- function(hexGrd = hex1,
                          n = 12,
                          subDat = DAT1) {
  rn <- 1:nrow(hexGrd)
  coords <- st_coordinates(st_centroid(hexGrd))
  k1 <- knn2nb(knearneigh(coords, k = n))
  all.linked <- max(unlist(nbdists(k1, coords)))
  nb.0.all <- dnearneigh(coords, 0, all.linked, row.names = rn)
  W.list <- nb2listw(nb.0.all, style = "B", zero.policy = TRUE)
  
  # Subset to remove missing cells
  id <- as.numeric(hexGrd$id)
  W.list <- subset(W.list,  id %in% unique(subDat$id))
  W <- nb2mat(nb.0.all, style = "B")
  W <- W[id %in% unique(subDat$id), id %in% unique(subDat$id)]
  row.names(W) <- id[id %in% unique(subDat$id)]
  
  # Plot to check
  fcell <- sample(1:length(nb.0.all), 1)
  hexGrd$id2 <- NA
  hexGrd$id2[nb.0.all[[fcell]]] <- 2
  hexGrd$id2[fcell] <- 1
  plot(hexGrd[, 10])
  
  return(list(W.list, W))
  
}

#' Check models (SAC)
modCheck <- function(modObj, maxCov) {
  
  #- Process data to grain and threshold
  modDat <-
    left_join(
      modObj$data[, c('covComp',
                      'area',
                      'accessibility',
                      'natureSpace',
                      'remoteness',
                      'nRecLog10',
                      'sppGrp')],
      maxCov,
      by = c('covComp',
             'area',
             'accessibility',
             'natureSpace',
             'remoteness',
             'nRecLog10',
             'sppGrp')
    )
  
  #- Check for spatial autocorrelation and model fit
  mod_pp <- posterior_predict(modObj)
  mod_ppm <- apply(mod_pp, MARGIN = 2, FUN = median)
  mod_res <- createDHARMa(
    simulatedResponse = t(mod_pp),
    fittedPredictedResponse = mod_ppm,
    observedResponse = modDat$covComp,
    integerResponse = TRUE
  )
  coord <- modDat[, c('X', 'Y', 'id')] %>% distinct()
  sa <-
    testSpatialAutocorrelation(
      recalculateResiduals(mod_res, 
                           group = modDat$id),
      x = coord[, 'X'],
      y = coord[, 'Y'])
  sa
  
}

#' Run bernoulli brms per species group
species.level.models <-
  function(hexCellArea,
           hexGrd,
           covEnvDat,
           nk = 12,
           sppGrpNames = unique(sppNameTab[, 2])) {
    # Model data
    DAT <- prepareModelDat(hexGrd, hexCellArea, covEnvDat)
    
    #~# Fit spatial model with species-level random effect
    for (sppGrp in sppGrpNames) {
      if (!file.exists(
        paste0(
          "Outputs/Figures_R1/brms_models/brms_",
          hexCellArea,
          "_CAR_",
          sppGrp,
          ".Rdata"
        )
      )) {
        # Subset species data
        sppDat <- DAT[DAT$sppGrp == sppGrp,]
        # Create neighbourhood matrix 3 deep nk = 12
        W <- create.nb.mat(hexGrd, n = nk, sppDat)
        # Fit model
        sppGrpMod <-
          brm(
            covComp |
              weights(area) ~ remoteness + accessibility +
              natureSpace + nRecLog10 + car(W, id, type = 'escar'),
            data = sppDat,
            data2 = list(W = W[[2]]),
            family = bernoulli(),
            chain = 2,
            iter = 5000,
            thin = 1,
            cores = 2,
            prior = set_prior(horseshoe(
              df = 3, par_ratio = 0.1
            )),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
          )
        # Check SA
        # moran.mc(residuals(sppGrpMod)[,1], listw = W[[1]], nsim = 99)
        # modCheck(sppGrpMod, sppDat)
        #set_prior(horseshoe(df = 3, par_ratio = 0.1))
        
        save(
          sppGrpMod,
          file = paste0(
            "Outputs/Figures_R1/brms_models/brms_",
            hexCellArea,
            "_CAR_",
            sppGrp,
            ".Rdata"
          ),
          compress = 'xz'
        )
        
      }
      
    }
    
    return(NULL)
    
  }

#' Extract species-group level model results
sppModResExtr <- function(sppGrp , grains) {
  # All model combinations
  sppMods <- expand.grid(sppGrp = sppGrp, grains = grains)
  # Loop through all combinations i <- 1
  resOutAll <- lapply(1:nrow(sppMods), function(i) {
    #- Load data
    modRes <-
      get(load(
        paste0(
          "Outputs/Figures_R1/brms_models/brms_",
          sppMods[i, 2],
          "_CAR_",
          sppMods[i, 1],
          ".Rdata"
        )
      ))
    #- Get R2 and summary
    R2 <- bayes_R2(modRes)
    modSum <- summary(modRes)
    #- Test specific hypothesis
    Remoteness_H1 <- hypothesis(modRes, 'remoteness < 0')
    NatureSpace_H1 <- hypothesis(modRes, 'natureSpace > 0')
    Accessibility_H1 <- hypothesis(modRes, 'accessibility > 0')
    NumberOfRecorders_H1 <- hypothesis(modRes, 'nRecLog10 > 0')
    
    #- Create output data.frame
    resOut <- data.frame(
      sppGrp = sppMods[i, 1],
      grain = sppMods[i, 2],
      nObs = modSum$nobs,
      # Intercept
      Intercept_b = modSum$fixed[1, 1],
      Intercept_L = modSum$fixed[1, 3],
      Intercept_U = modSum$fixed[1, 4],
      # Remoteness
      Remoteness_b = modSum$fixed[2, 1],
      Remoteness_L = modSum$fixed[2, 3],
      Remoteness_U = modSum$fixed[2, 4],
      Remoteness_H1_ER = Remoteness_H1$hypothesis$Evid.Ratio,
      Remoteness_H1_PProb = Remoteness_H1$hypothesis$Post.Prob,
      # Access
      Accessibility_b = modSum$fixed[3, 1],
      Accessibility_L = modSum$fixed[3, 3],
      Accessibility_U = modSum$fixed[3, 4],
      Accessibility_H1_ER = Accessibility_H1$hypothesis$Evid.Ratio,
      Accessibility_H1_PProb = Accessibility_H1$hypothesis$Post.Prob,
      # Green space
      NatureSpace_b = modSum$fixed[4, 1],
      NatureSpace_L = modSum$fixed[4, 3],
      NatureSpace_U = modSum$fixed[4, 4],
      NatureSpace_H1_ER = NatureSpace_H1$hypothesis$Evid.Ratio,
      NatureSpace_H1_PProb = NatureSpace_H1$hypothesis$Post.Prob,
      # N recorders
      NumberOfRecorders_b = modSum$fixed[5, 1],
      NumberOfRecorders_L = modSum$fixed[5, 3],
      NumberOfRecorders_U = modSum$fixed[5, 4],
      NumberOfRecorders_H1_ER = NumberOfRecorders_H1$hypothesis$Evid.Ratio,
      NumberOfRecorders_H1_PProb = NumberOfRecorders_H1$hypothesis$Post.Prob,
      # R2
      R2_b = R2[1, 1],
      R2_L = R2[1, 3],
      R2_U = R2[1, 4]
    )
    
    return(resOut)
    
  })
  resOutAll <- do.call(bind_rows, resOutAll)
  return(resOutAll)
  
}

#' Run zeta distance decay over time
temporal.zeta <- function(bioDat,
                          sppGrpName,
                          timeSlice = 5,
                          region_sf = cornwall_sf,
                          grain = grn,
                          hexGrid = hexGrain) {
  #~# Select all years
  siteXSppYr <-
    siteXSpeciesXTime(
      bioDat,
      timeSlice,
      region_sf,
      hexCellArea = grain,
      bioHex = NULL,
      hexGrid = hexGrid
    )
  siteXSppYr$site <-
    as.numeric(sapply(siteXSppYr$siteTime, function(x)
      strsplit(x, split = "_")[[1]][1]))
  siteXSppYr$year <-
    sapply(siteXSppYr$siteTime, function(x)
      strsplit(x, split = "_")[[1]][2])
  
  #~# Run zeta for timeslices # si <- 8
  sites <- mixedsort(unique(siteXSppYr$site))
  
  sfInit(parallel = TRUE, cpus = 7)
  sfExport(list = c("siteXSppYr"))
  sfLibrary(dplyr)
  sfLibrary(zetadiv)
  
  sampCov <- sfLapply(sites, function(si) {
    print(si)
    #~# Sample increasing amounts of time
    ts_dat_tmp <- siteXSppYr %>%
      filter(site == si) %>%
      select(!c('siteTime', 'site')) %>%
      group_by(year) %>%
      summarise(across(.fns = max), .groups = 'drop') %>%
      filter(year >= 1975)
    #~# xy
    xy <-
      data.frame(x = rep(1, 9),
                 y = 1:9,
                 year = seq(1975, 2015, 5))
    xy <- xy[xy$year %in% ts_dat_tmp$year, 1:2]
    #~# Remove year
    ts_dat_tmp <-
      ts_dat_tmp %>% select(!c('year')) %>% as.data.frame()
    #~# Run zeta decline
    possibleError <- tryCatch(
      zeta_dd <-
        Zeta.ddecay(
          xy = xy,
          data = ts_dat_tmp,
          order = 2,
          reg.type = 'glm',
          normalize = 'Simpson',
          plot = F
        ),
      error = function(e)
        e
    )
    
    if (!inherits(possibleError, "error")) {
      zeta_dd <- summary(zeta_dd$reg)
      dd_est <-  round(zeta_dd$coefficients[2, 1], 3)
      dd_se <-  round(zeta_dd$coefficients[2, 2], 3)
      dd_t <-  round(zeta_dd$coefficients[2, 3], 3)
      dd_p <-  round(zeta_dd$coefficients[2, 4], 3)
      resOut <-
        data.frame(
          site = si,
          dd_est = dd_est,
          dd_se = dd_se,
          dd_t = dd_t,
          dd_p = dd_p
        )
    } else {
      resOut <-
        data.frame(
          site = si,
          dd_est = NA,
          dd_se = NA,
          dd_t = NA,
          dd_p = NA
        )
    }
    return(resOut)
    
  })
  sampCov <- do.call(rbind, sampCov)
  sampCov$sppGrp <- sppGrpName
  sampCov$grain <- grain
  sampCov$nCell <- length(unique(siteXSppYr$site))
  sampCov
  
}


#' Simpson's diversity index
simpsonIndex <- function(x) {
  N <- sum(x)
  p <- (x / N) ^ 2
  D <- 1 / sum(p)
  
}

