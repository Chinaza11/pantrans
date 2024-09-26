# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 06:10:21 2024

@author: Chinaza
"""

import cv2, os

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# Read the two images
image1 = cv2.imread('Rio_reference_genome/orthofinder_3/rio_2.png')
image2 = cv2.imread('pantranscriptome/orthofinder_3/pant_2.png')

# Resize image2 to have the same width as image1 (for vertical concatenation)
if image1.shape[1] != image2.shape[1]:
    image2 = cv2.resize(image2, (image1.shape[1], image2.shape[0]))

# Vertical concatenation: place image1 above image2
combined_vertical = cv2.vconcat([image1, image2])

# Save the combined image
cv2.imwrite('final_processing_and_plotting/duplication_event_enrichment_plot.png', combined_vertical)