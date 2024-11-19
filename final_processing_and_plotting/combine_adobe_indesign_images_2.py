# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 22:31:07 2024

@author: Chinaza
"""

import os
from PIL import Image

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

def combine_images(image1_path, image2_path, combined_image_path, dpi=300, width_cm=20):
    # Open the two images using PIL
    image1 = Image.open(image1_path)
    image2 = Image.open(image2_path)

    # Resize image2 to have the same width as image1 (for vertical concatenation)
    if image1.width != image2.width:
        image2 = image2.resize((image1.width, image2.height))

    # Create a new image with the height of both images combined
    total_height = image1.height + image2.height
    combined_image = Image.new('RGB', (image1.width, total_height))

    # Paste image1 and image2 into the combined image
    combined_image.paste(image1, (0, 0))  # Paste image1 at the top
    combined_image.paste(image2, (0, image1.height))  # Paste image2 below image1

    # Convert width from cm to pixels
    width_px = int((width_cm * dpi) / 2.54)

    # Calculate the new height to maintain aspect ratio
    aspect_ratio = combined_image.height / combined_image.width
    new_height_px = int(width_px * aspect_ratio)

    # Resize the concatenated image to the desired width (20 cm) while maintaining aspect ratio
    resized_combined_image = combined_image.resize((width_px, new_height_px))

    # Save the combined image as a PDF with the specified DPI
    resized_combined_image.save(combined_image_path, "PDF", resolution=dpi)

# =============================================================================
# Main result
# =============================================================================

image1_path = 'Rio_reference_genome/orthofinder_3/rio_2_2.jpg'
image2_path = 'pantranscriptome/orthofinder_3/pant_2_2.jpg'

combined_image_path = 'final_processing_and_plotting/duplication_event_enrichment_plot_2.pdf'

combine_images(image1_path, image2_path, combined_image_path)

# =============================================================================
# permutation result
# =============================================================================

image1_path = 'Rio_reference_genome/orthofinder_3/permutation/rio_permutation.png'
image2_path = 'pantranscriptome/orthofinder_3/permutation/pant_permutation.png'
combined_image_path = 'final_processing_and_plotting/duplication_event_enrichment_plot-permutation_result_2.pdf'

combine_images(image1_path, image2_path, combined_image_path)