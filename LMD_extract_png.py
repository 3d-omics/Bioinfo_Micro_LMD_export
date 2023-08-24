# Extract PNG images from an RTF file
# By: Antton Alberdi (antton.alberdi@gmail.com)
# Date: 12/08/2023
# Usage: python extract_png.py -r input/D015cH109_230811.rtf
# Output: png images inside the "images" folder

import re
import binascii
import argparse

parser = argparse.ArgumentParser(description="Extract PNG images from an RTF file.")
parser.add_argument("-r", "--rtf_file", required=True, help="Path to the RTF file containing the image data")
args = parser.parse_args()

# Read the RTF file containing the image data
rtf_file = args.rtf_file
with open(rtf_file, "r") as f:
    rtf_text = f.read()

# Find all occurrences of PNG binary data using regex
png_matches = re.findall(r"\\pngblip(.*?)\\par", rtf_text, re.DOTALL)

# Find image names using regex
image_names = re.findall(r"\\f2 {\\ltrch (.*?)}\\li0\\ri0\\sa0\\sb0\\fi0\\qj\\par}\n{\\f2 {\\ltrch {\\\*\\shppict{\\pict", rtf_text)

# Process each PNG match and corresponding image name
for idx, (png_match, image_name) in enumerate(zip(png_matches, image_names), start=1):
    png_data_hex = re.sub(r"[^0-9a-fA-F]", "", png_match)  # Clean non-hex characters
    png_data = binascii.unhexlify(png_data_hex)
    output_png_file = "images/{}.png".format(image_name)
    with open(output_png_file, "wb") as f:
        f.write(png_data)
    print("PNG image {} saved as {}".format(idx, output_png_file))
