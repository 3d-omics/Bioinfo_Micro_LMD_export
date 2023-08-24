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
git clone https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics.git
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

Rscript LMD_report.R \
 -b {batchname}
 -r {trf input file(s)}
 -c {csv input file}
 -t {tsv input file}
 -p {pdf input file}
 -w {html input file}

If using multiple RTF files, list them comma-separated without spaces as shown below.
```{sh}
Rscript LMD_report.R -b D015cH109_230811 -r input/D015cH109_230811_1.rtf,input/D015cH109_230811_2.rtf,input/D015cH109_230811_3.rtf -c input/D015cH109_230811.csv -t output/D015cH109_230811.tsv -p output/D015cH109_230811.pdf -w output/D015cH109_230811.html
```

### 4. Check the output
The following output files should be found in the output folder:
- output/D015cH109_230811.html
- output/D015cH109_230811.pdf
- output/D015cH109_230811.tsv
