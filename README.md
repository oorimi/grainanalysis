# grainanalysis
Matlab script to analyze the grain size from a grain boundary image

## Introduction
This repository provides two MATLAB scripts: "imgseg.m" that contains functions for image processing, and "demo.m" that runs the analysis with user-specified parameters. Given an image where grain boundaries are traced, the scripts extract the size of each grain. A grain boundary pixel is merged randomly into one of the neighboring grains if it is shared by more than one grain. The colors in the grain map are randomly generated. One example is given in "example_data" and the corresponding results can be found in "example_results".

## Steps to use:
1. Ensure the image to process (e.g. "example.jpg") and the codes are under the same directory
2. Open "demo.m" in MATLAB and edit the following information:
   ```
       imageName: name of the image
       threshold: threshold value for binarizing the image
       pixel2um: scale bar of the image, i.e. number of pixels per micronmeter
       savename: name for saving the results   
   ```
3. Run "demo.m"
4. Check the results summary printed in the MATLAB command line and the files saved in the same directory
