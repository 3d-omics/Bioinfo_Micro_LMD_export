#' ---
#' title: "LMD report"
#' author: "`r params$batch`"
#' output:
#'   pdf_document:
#'     toc: true
#' params:
#'   batch: NULL
#'   rtf: NULL
#'   csv: NULL
#'   tsv: NULL
#' ---

#'
#' # Introduction
#' This is the laser microdissection report automatically generated based on the the LMD7 output files.
#' Here you will find basic information about the microdissected samples.
#' The reporting script is still in developmental phase.

#+ load_libraries, echo=FALSE, warning=FALSE
#Load libraries
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(striprtf))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(png))

#+ slide_positions, echo=FALSE, warning=FALSE
#Based on Bryan measurements in 2023/08/11
#tl=top-left, br=bottom-right
slide1_tl <- c(250,2212) # (x, y) for top-left
slide1_br <- c(25660,77470)    # (x, y) for bottom-right
membrane1_tl <- c(4962,8698) # (x, y) for top-left
membrane1_br <- c(21013,53732) # (x, y) for bottom-right

slide2_tl <- c(28983,2004) # (x, y) for top-left
slide2_br <- c(54565,77541) # (x, y) for bottom-right
membrane2_tl <- c(33689,8698) # (x, y) for top-left
membrane2_br <- c(49650,53658) # (x, y) for bottom-right

slide3_tl <- c(64978,1938) # (x, y) for top-left
slide3_br <- c(90644,77372) # (x, y) for bottom-right
membrane3_tl <- c(69768,8658) # (x, y) for top-left
membrane3_br <- c(85819,53692) # (x, y) for bottom-right

slide4_tl <- c(93780,1946) # (x, y) for top-left
slide4_br <- c(119690,77465) # (x, y) for bottom-right
membrane4_tl <- c(98718,8682) # (x, y) for top-left
membrane4_br <- c(114760,53772) # (x, y) for bottom-right

#+ load_files, echo=FALSE, warning=FALSE
# Load input files

#batch<-"230823RTest_1"
#csv_file="input/230823RTest1.csv"
#rtf_files="input/230823RTest_1.rtf,input/230823RTest_2.rtf,input/230823RTest_3.rtf,input/230823RTest_4.rtf"

batch <- params$batch
csv <- read.csv(params$csv, sep=";")

#Load from single or multiple RTF files
if (length(strsplit(params$rtf, ",")[[1]]) == 1) {
  rtf <- read_rtf(params$rtf)
} else {
  rtf <- strsplit(params$rtf, ",")[[1]] %>%
      lapply(., read_rtf) %>%
      unlist()
}

#+ rtf_extract_text, echo=FALSE, warning=FALSE
#Field of view image name
imageID <- grep("LMDPic", rtf, value = TRUE)

#Field of view image width
imageDimX <- grep("[0-9][0-9]\\s\\(xy\\)", rtf, value = TRUE) %>%
    str_match(., "(\\d+) x \\d+ \\(xy\\)") %>%
    as.data.frame() %>%
    select(2) %>%
    pull() %>%
    as.numeric()
#Field of view image height
imageDimY <- grep("[0-9][0-9]\\s\\(xy\\)", rtf, value = TRUE) %>%
    str_match(., "x (\\d+) \\(xy\\)") %>%
    as.data.frame() %>%
    select(2) %>%
    pull() %>%
    as.numeric()
#Grid size (dimension per pixel)
grid <- grep("Voxel Size:", rtf, value = TRUE) %>%
  str_match(., "Size: ([^ ]+)") %>%
  as.data.frame() %>%
  select(2) %>%
  pull() %>%
  as.numeric()
#Grid unit (nm, um, etc.)
gridunit <- grep("Voxel Size:", rtf, value = TRUE) %>%
  str_match(., "\\S+(?= \\(xy\\))") %>%
  c()
#Convert all grids to µm
conversion_factors <- c("nm" = 0.001, "µm" = 1, "mm" = 1000)
desired_magnitude <- "µm"
grid <- grid * conversion_factors[gridunit] / conversion_factors[desired_magnitude]
#Date and time of the image
magnification <- grep("Magnification", rtf, value = TRUE) %>%
  gsub("Magnification:","",.)
#Date and time of the image
timestamp <- grep("Date", rtf, value = TRUE) %>%
  gsub("Date: ","",.)
#X axis of stage position in um
stagePositionX <- grep("Stage", rtf, value = TRUE) %>%
  str_match(., "X:(\\d+)") %>%
  as.data.frame() %>%
  select(2) %>%
  pull() %>%
  as.numeric() / 10 #Convert in to microns
#Y axis of stage position in um
stagePositionY <- grep("Stage", rtf, value = TRUE) %>%
  str_match(., "Y:(\\d+)") %>%
  as.data.frame() %>%
  select(2) %>%
  pull() %>%
  as.numeric() / 10 #Convert in to microns
#Compile all information into a table
LMDPic_table <- data.frame(imageID,imageDimX,imageDimY,timestamp,stagePositionX,stagePositionY,grid,magnification) %>%
  filter(!grepl("_ins_", imageID)) #remove collector images

#+ rtf_extract_shapes, echo=FALSE, warning=FALSE
# Grep shape lines
shapes <- grep("Shape", rtf, value = TRUE)
# Replace "ShapeData:" by imageIDs
shapedata_indices <- which(shapes == "ShapeData:")
for (i in seq_along(shapedata_indices)) {
  shapes[shapedata_indices[i]] <- LMDPic_table$imageID[i]
}
# Fetch Collector and Area information from each shape
original_vector <- shapes
current_string <- NULL
matrix_rows <- list()
for (element in original_vector) {
  if (grepl("LMDPic", element)) {
    current_string <- element
  } else if (!is.null(current_string)) {
    collector <- str_match(element, "Collector: ([^ ]+)")[,2]
    area <- str_match(element, "Area: (\\d+)")[,2]
    matrix_rows <- append(matrix_rows, list(c(current_string, collector, area)))
  }
}
# Create shape table
shape_table <- do.call(rbind, matrix_rows) %>%
  as.data.frame() %>%
  rename(imageID=V1,collector=V2,area=V3) %>%
  filter(grepl("_pre_", imageID)) %>% #only retain pre-cut images
  filter(area > 0) #only retain rows with cut areas >0 (0s are usually lines used for corrections)

#+ detect_collector_type, echo=FALSE, warning=FALSE
if(nchar(shape_table[1,2]) == 1){
  collector_type="strip"
}else{
  collector_type="plate"
  }

#+ rtf_merge_shapes, echo=FALSE, warning=FALSE
if(collector_type == "plate"){
  rtf_table <- inner_join(LMDPic_table,shape_table,by="imageID") %>%
    separate(collector, into = c("row", "col"), sep = "(?<=\\D)(?=\\d)", remove = FALSE) %>%
    mutate(row = ifelse(row == "ParkPosition", NA, row))
}

if(collector_type == "strip"){
  rtf_table <- inner_join(LMDPic_table,shape_table,by="imageID") %>%
    rename(row=collector) %>%
    mutate(area = ifelse(is.na(as.numeric(area)), 0, as.numeric(area))) %>%
    mutate(col = ifelse(area != lag(area, default = -1), cumsum(area != lag(area, default = -2)), cumsum(area != lag(area, default = -1)))) %>% #infer column
    mutate(collector = paste(row, col, sep = "")) %>%
    mutate(row = ifelse(row == "ParkPosition", NA, row))
}

#+ csv_transform, echo=FALSE, warning=FALSE
#Define column name extraction function
extract_colname <- function(...) {
  colnames <- c("A", "B", "C", "D", "E", "F", "G", "H")
  non_na_col <- colnames[!is.na(c(...))]
  if (length(non_na_col) > 0) {
    return(paste(non_na_col, collapse = ", "))
  } else {
    return("")
  }
}
#Transform csv
csv_table <- csv %>%
  rename(A=A..Area., B=B..Area., C=C..Area., D=D..Area., E=E..Area., F=F..Area., G=G..Area., H=H..Area., type=Type, imageID2=Image.Name.s., coord=X.Y.Coordinates) %>%
  slice(-1) %>%
  separate(coord, into = c("pixelX", "pixelY"), sep = " / ") %>% #split x and y coordinates
  mutate(pixelX = ifelse(is.na(as.numeric(pixelX)), 0, as.numeric(pixelX))) %>%
  mutate(pixelY = ifelse(is.na(as.numeric(pixelY)), 0, as.numeric(pixelY))) %>%
  mutate(area = rowSums(.[, c("A", "B", "C", "D", "E", "F", "G", "H")], na.rm = TRUE)) %>% #transfer area information to new column
  mutate(row = purrr::pmap_chr(select(., A:H), extract_colname)) %>% #infer row letter from position
  select(imageID2, type, pixelX, pixelY, area, row) %>% #subset relevant columns %>%
  filter(type == "Ellipse") %>%
  mutate(col = ifelse(area != lag(area, default = -1), cumsum(area != lag(area, default = -2)), cumsum(area != lag(area, default = -1)))) %>% #infer column
  mutate(collector = paste(row, col, sep = ""))

#+ rtf_csv_merge, echo=FALSE, warning=FALSE
full_table <- rtf_table %>%
  filter(area != 0) %>%
  filter(collector != "ParkPosition") %>%
  select(stagePositionX,stagePositionY,imageID,imageDimX,imageDimY,collector,grid,magnification) %>%
  group_by(collector) %>%
  slice(n()) %>%
  ungroup() %>%
  left_join(csv_table, by = join_by(collector==collector)) %>%
  arrange(col,row) %>%
  filter(!duplicated(paste(collector, area, sep = "_"), fromLast = TRUE)) #select the last appearance whenver there are repetitions

#+ calculate_positions, echo=FALSE, warning=FALSE
position_table <- full_table %>%
  mutate(stagePositionX = ifelse(is.na(as.numeric(stagePositionX)), 0, as.numeric(stagePositionX))) %>%
  mutate(stagePositionY = ifelse(is.na(as.numeric(stagePositionY)), 0, as.numeric(stagePositionY))) %>%
  mutate(imageDimX = ifelse(is.na(as.numeric(imageDimX)), 0, as.numeric(imageDimX))) %>%
  mutate(imageDimY = ifelse(is.na(as.numeric(imageDimY)), 0, as.numeric(imageDimY))) %>%
  mutate(pixelX = ifelse(is.na(as.numeric(pixelX)), 0, as.numeric(pixelX))) %>%
  mutate(pixelY = ifelse(is.na(as.numeric(pixelY)), 0, as.numeric(pixelY))) %>%
  mutate(grid = ifelse(is.na(as.numeric(grid)), 0, as.numeric(grid))) %>%
  # Calculate exact cut positions (center of ellipse)
  mutate(cutPositionX=stagePositionX-(imageDimX*grid/2)+pixelX*grid) %>%
  mutate(cutPositionY=stagePositionY-(imageDimY*grid/2)+pixelY*grid) %>%
  # Calculate edge positions of the FOV images
  mutate(fieldX=stagePositionX-(imageDimX*grid/2)) %>%
  mutate(fieldY=stagePositionY-(imageDimY*grid/2)) %>%
  # Calculate width and height of FOV images
  mutate(fieldW=imageDimX*grid) %>%
  mutate(fieldH=imageDimY*grid) %>%
  mutate(area = ifelse(is.na(as.numeric(area)), 0, as.numeric(area))) %>%
  # Calculate radius of ellipse from area
  mutate(radius=round(sqrt(area / pi),2)) %>%
  select(collector, cutPositionX, cutPositionY, imageID, fieldX, fieldY, fieldW, fieldH, radius, col, area, magnification)

#+ position_table_print, echo=FALSE, warning=FALSE
position_table %>%
  select(collector, cutPositionX, cutPositionY,radius, area, fieldX, fieldY, fieldW, fieldH) %>%
  kable()

#+ position_table_output, echo=FALSE, warning=FALSE
position_table %>%
  select(collector, cutPositionX, cutPositionY,radius, area, fieldX, fieldY, fieldW, fieldH) %>%
  write.table(., file=params$tsv, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


#+ extract_png, echo=FALSE, warning=FALSE, comments="", message=FALSE, results="hide"
# Define python script
launch_python_script <- function(value) {
  system(paste0("python LMD_extract_png.py -r ", value), wait=TRUE)
}

#Launch python script
if (length(strsplit(params$rtf, ",")[[1]]) == 1) {
  launch_python_script(params$rtf)
} else {
  strsplit(params$rtf, ",")[[1]] %>%
      lapply(., launch_python_script)
}

#+ load_png, echo=FALSE, warning=FALSE
png_files <- list.files(path = "images", pattern = "\\.png$", full.names = TRUE)
png_names <- tools::file_path_sans_ext(basename(png_files))
png_img <- map(png_files, readPNG)
png_tibble <- tibble(imageID = png_names, image = png_img)
#Merge to position_table
position_img_table <- position_table %>%
  left_join(png_tibble, by = join_by(imageID==imageID))

#+ cutting_overview, echo=FALSE, warning=FALSE
# Plot height
height=80000
# Slide background
overview_plot <- ggplot() +
  #draw slides
  geom_rect(aes(xmin = slide1_tl[1], xmax = slide1_br[1], ymin = slide1_tl[2], ymax = slide1_br[2]), fill = "lightgrey", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = slide2_tl[1], xmax = slide2_br[1], ymin = slide2_tl[2], ymax = slide2_br[2]), fill = "lightgrey", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = slide3_tl[1], xmax = slide3_br[1], ymin = slide3_tl[2], ymax = slide3_br[2]), fill = "lightgrey", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = slide4_tl[1], xmax = slide4_br[1], ymin = slide4_tl[2], ymax = slide4_br[2]), fill = "lightgrey", color = "grey", lwd=0.1) +
  #draw membranes
  geom_rect(aes(xmin = membrane1_tl[1], xmax = membrane1_br[1], ymin = membrane1_tl[2], ymax = membrane1_br[2]), fill = "#f4f4f4", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = membrane2_tl[1], xmax = membrane2_br[1], ymin = membrane2_tl[2], ymax = membrane2_br[2]), fill = "#f4f4f4", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = membrane3_tl[1], xmax = membrane3_br[1], ymin = membrane3_tl[2], ymax = membrane3_br[2]), fill = "#f4f4f4", color = "grey", lwd=0.1) +
  geom_rect(aes(xmin = membrane4_tl[1], xmax = membrane4_br[1], ymin = membrane4_tl[2], ymax = membrane4_br[2]), fill = "#f4f4f4", color = "grey", lwd=0.1)


    # Add overview rasters
    for (i in 1:nrow(position_img_table)){
      overview_plot <- overview_plot +
        annotation_raster(position_img_table$image[[i]], xmin = position_img_table$fieldX[[i]], xmax = position_img_table$fieldX[[i]] + position_img_table$fieldW[[i]], ymax = -position_img_table$fieldY[[i]] - position_img_table$fieldH[[i]], ymin = -position_img_table$fieldY[[i]])
    }

    # Add cutting ellipses and style
    overview_plot <-  overview_plot +
      #draw views (for some reason, erases the pictures)
      #geom_rect(data = position_table, aes(xmin = fieldX, xmax = fieldX+fieldW, ymin = height-fieldY-fieldH, ymax = height-fieldY, color = as.factor(col)), fill = NA, lwd=0.1)+
      #draw cuts
      geom_circle(data = position_table, aes(x0=cutPositionX, y0=cutPositionY, r = radius, fill=as.factor(col), color=as.factor(col), linetype="blank"),alpha=0.8) +
      coord_fixed() +
      xlim(0, 120000) +
      ylim(80000, 0) +
      labs(title = "Cut overview", x = NULL, y = NULL) +
      theme_classic()

  #+ cutting_overview_plot, echo=FALSE, warning=FALSE, fig.height=7
  overview_plot

  #' \newpage

  #+ cutting_single_plot, echo=FALSE, warning=FALSE, fig.height=8
  for (i in 1:nrow(position_img_table)){
    #Subset table
    single <- position_img_table %>% slice(i)
    #Generate image
    single_plot <- ggplot() +
      xlim(single$fieldX, single$fieldX + single$fieldW) +
      ylim(single$fieldY + single$fieldH, single$fieldY) +
      coord_fixed() +
      annotation_raster(single$image[[1]], xmin = single$fieldX, xmax = single$fieldX + single$fieldW, ymin = -single$fieldY, ymax = -single$fieldY - single$fieldH) +
      geom_circle(data = single, aes(x0=cutPositionX, y0=cutPositionY, r = radius, fill=as.factor(col), color=as.factor(col), linetype="blank"),alpha=0.8) +
      geom_vline(data = single, aes(xintercept = cutPositionX, color=as.factor(col)), alpha=0.5)+
      geom_hline(data = single, aes(yintercept = cutPositionY, color=as.factor(col)), alpha=0.5)+
      theme_classic() +
      theme(legend.position="none", plot.margin = margin(t = 7, unit = "lines"))

    #Create above text
    single_collector <- paste0("Collector: ",single$collector)
    single_coords <- paste0("X-axis: ",single$cutPositionX," µm\nY-axis: ",single$cutPositionY," µm\nArea: ",single$area," µm\nDiameter: ",single$radius *2," µm\nMagnification: ",single$magnification,"X")

    #Merge text with figure
    single_plot_annot <-  annotate_figure(single_plot,
          top = arrangeGrob(
              text_grob(single_collector, , x = unit(0, "npc"), y = unit(1, "npc") - unit(0.5, "lines"), size = 14, face = "bold", just = "left"),
              text_grob(single_coords, x = unit(0, "npc"), y = unit(1, "npc") - unit(4, "lines"), size = 12, just = "left"),
              ncol = 1)
      )
    #Print merged image
    print(single_plot_annot)
  }
