####
# LMD report rendering script
# By: Antton Alberdi (antton.alberdi@sund.ku.dk)
# Date: 24/08/2023
# Dependencies: pandoc
####

#Load libraries
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(R.utils))

#Parse input-output
parser <- ArgumentParser()
parser$add_argument("-b", "--batch", action="store", help="Batch/experiment code")
parser$add_argument("-r", "--rtf", action="store", help="RTF file")
parser$add_argument("-c", "--csv", action="store", help="CSV file")
parser$add_argument("-t", "--tsv", action="store", help="TSV output file")
parser$add_argument("-p", "--pdf", action="store", help="PDF output file")
parser$add_argument("-w", "--html", action="store", help="HTML output file")
args <- parser$parse_args()

#Create required directories
dir.create("images")
dir.create("output")

#Render pdf document using LMD_report_v2
rmarkdown::render("LMD_report_v2.R", params = list(batch=args$batch, rtf = args$rtf, csv = args$csv, tsv = args$tsv), output_format = "pdf_document", output_file = args$pdf)

#Render html document using LMD_report_v2
rmarkdown::render("LMD_report_v2.R", params = list(batch=args$batch, rtf = args$rtf, csv =args$csv, tsv = args$tsv), output_format = "html_document", output_file = args$html)

#cd /Users/anttonalberdi/github/Bioinfo_Micro_LMD_data_export/
#Rscript LMD_report.R -b D015cH109 -r input/D015cH109_230811.rtf -c input/D015cH109_230811.csv -p output/D015cH109_230811.pdf -w output/D015cH109_230811.html
#Rscript LMD_report.R -b 230823RTest_1 -r input/230823RTest_1.rtf,input/230823RTest_2.rtf,input/230823RTest_3.rtf,input/230823RTest_4.rtf -c input/230823RTest1.csv -t output/230823RTest1.tsv -p output/230823RTest1.pdf -w output/230823RTest1.html
#Rscript LMD_report.R -b 230823RTest_2 -r input/230823RTest2_1.rtf,input/230823RTest2_2.rtf,input/230823RTest2_3.rtf -c input/230823RTest2.csv -t output/230823RTest1.tsv -p output/230823RTest2.pdf -w output/230823RTest2.html
