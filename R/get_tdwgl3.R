#' @title Get TDWG World Geographical Scheme for Recording Plant Distributions (WGSRPD) data
#'
#' @description Function to get the WGSRPD (TDWG) spatial data
#' @keywords TDWG WGSRPD
#' @import sf
#' @return geojson object with TDWG level 3 boundaries
#' @export


# changed this to shapefile as geojson is too coarse
# problem with downloads - doesn't create valid shp.


get_tdwgl3 <- function() {
  if (!file.exists("level3.shp")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.shp",
      "data-raw/level3.shp",
          )
   }
  if (!file.exists("level3.dbf")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.dbf",
      "data-raw/level3.dbf"
    )
    
  }
  if (!file.exists("level3.sbn")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.sbn",
      "data-raw/level3.sbn"
    )
    
  }
  if (!file.exists("level3.sbx")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.sbx",
      "data-raw/level3.sbx"
    )
    
  }
  if (!file.exists("level3.prj")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.prj",
      "data-raw/level3.prj"
    )
    
  }
  if (!file.exists("level3.shx")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.shx",
      "data-raw/level3.shx"
    )
    
  }
  if (!file.exists("level3.shp.xml")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.prj",
      "data-raw/level3.shp.xml"
    )
    
  }
  
  if (!file.exists("level3.level4_Lev.atx")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/level3.level4_Lev.atx",
      "data-raw/level3.level4_Lev.atx"
    )
    
  }
  if (!file.exists("master_meta_data.html")) {
    download.file(
      "https://github.com/tdwg/wgsrpd/raw/master/level3/master_meta_data.html",
      "data-raw/master_meta_data.html"
    )
    
  }
  
  
}

