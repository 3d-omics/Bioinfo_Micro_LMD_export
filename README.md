# Bioinfo_Micro_LMD_data_export
This repository contains R and Python scripts for generating precise spacial microdissection information from raw CSV and RTF file produced by the Leica LMD7 software.

## Dependencies
- R
  - knitr
  - striprtf
  - tidyverse
  - ggplot2
  - ggforce
  - ggpubr
  - gridExtra
  - png
- Python
  - re
  - binascii
  - argparse

## How to use it

### 1. Prepare working directory
```{sh}
# Clone the repository
git clone https://github.com/3d-omics/Bioinfo_Micro_LMD_export.git
cd Bioinfo_Macro_Host_Transcriptomics
mkdir input #for storing the input RTF and CSV files
mkdir output #for storing output PDF, HTML and TSV files
mkdir images #for storing PNG images extracted from the RTF files
```
### 2. Add input data
Place desired RTF (could be multiple) and CSV (must be a single one) files in the input folder. They should look something like the following:
- input/D015cH109_230811_1.rtf
- input/D015cH109_230811_2.rtf
- input/D015cH109_230811_3.rtf
- input/D015cH109_230811.csv

### 3. Launch script

Rscript LMD_report.R
 - -b {batch name}
 - -r {rtf input file(s)}
 - -c {csv input file}
 - -t {tsv output file}
 - -p {pdf output file}
 - -w {html output file}

If using multiple RTF files, list them comma-separated without spaces as shown below.
```{sh}
Rscript LMD_report.R -b D015cH109_230811 -r input/D015cH109_230811_1.rtf,input/D015cH109_230811_2.rtf,input/D015cH109_230811_3.rtf -c input/D015cH109_230811.csv -t output/D015cH109_230811.tsv -p output/D015cH109_230811.pdf -w output/D015cH109_230811.html
```

### 4. Check the output
The following output files should be found in the output folder:
- output/D015cH109_230811.html
- output/D015cH109_230811.pdf
- output/D015cH109_230811.tsv

#### Positioning data

|collector|cutPositionX|cutPositionY|radius|area|fieldX|fieldY|fieldW|fieldH|
| ---  | ---  | ---  | ---  | ---  | ---  | ---  | ---  | ---  |
|A1 | 44627.28 | 13486.46 | 287.00 | 258769 | 42361.80 | 11852.60 | 5011.2000 | 3758.4000 |
|B1 | 44021.76 | 14324.27 | 287.00 | 258769 | 42361.80 | 11852.60 | 5011.2000 | 3758.4000 |
|C1 | 11893.07 | 21371.01 | 287.00 | 258769 | 9395.30 | 18993.30 | 5011.2000 | 3758.4000 |
|B2 | 75203.29 | 38863.49 | 19.79 | 1231 | 75124.60 | 38836.15 | 198.7968 | 149.0976 |
|C2 | 75221.83 | 38899.52 | 19.79 | 1231 | 75124.60 | 38836.15 | 198.7968 | 149.0976 |
|D2 | 75196.98 | 38950.46 | 19.79 | 1231 | 75124.60 | 38836.15 | 198.7968 | 149.0976 |
|E2 | 75268.00 | 38876.74 | 19.79 | 1231 | 75124.60 | 38836.15 | 198.7968 | 149.0976 |
|A3 | 11463.57 | 29818.68 | 176.10 | 97424 | 10443.42 | 29456.67 | 1252.3584 | 939.2688 |
|B3 | 45397.72 | 13727.23 | 176.10 | 97424 | 44867.42 | 13262.17 | 1252.3584 | 939.2688 |
|D3 | 11445.96 | 30196.34 | 176.10 | 97424 | 10443.42 | 29456.67 | 1252.3584 | 939.268 |
