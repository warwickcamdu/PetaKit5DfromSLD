[![DOI](https://zenodo.org/badge/829345148.svg)](https://zenodo.org/doi/10.5281/zenodo.12805990)

# PetaKit5DfromSLD

Modified version of PetaKit5D to read in the 3i Lattice Lightsheet Microscope Slidebook files (.sld), deconvolve and deskew using PetaKit5D, then reconstruct to a 5D tif. Requires installing PetaKit5D.

Authors: David Corcoran, Laura Cooper, Scott Brooks and Yara Aghabi

## Usage

Run ``` modularPipeline``` - Change to CAMDU_decon_wrapper

A folder selector window will open, please select a folder containing solely your PSFs

Your PSFs **MUST** be named in the following format:

PSF_CH0.tif 

PSF_CH1.tif

Up to a maximum of 9. Your first channel will be CH0.

We recommend that you include a readme file in this folder that documents any metadata lost in renaming the files i.e. Channel frequency/colour of PSFs.

Within the input folder you need to have a PSF file for each channel in .tif format (not .tiff). The PSFs must be imaged using "sample scan" and with the same Z-spacing as the images (e.g. 0.5um). The metadata needs to be correct for the XYZ pixel spacing and the units written as "um" not "microns" (e.g. 0.104 um for XY and 0.5 um for Z). The image must cropped to contain only one bead.

**TODO: Select PSF folder then allow user to pick which PSF corresponds to which channel order**

Once selected a second folder selector will open, please select the folder containing the input images.

All the .sld files to be processed should be placed in the same folder (not within subfolders). Each .sld file can contain multiple series (otherwise known as captures). Set the input folder to the folder containing the .sld files. Don't use "C0" or "C1" in the .sld filenames or anywhere in the pathname, otherwise it'll break. 




## Configuration

**TODO: Put common settings into a UI**


Processing mode: Choose 'deskew-only', 'decon+deskew', or 'both'.
```config.processingMode = 'both';```

Number of iterations for deconvolution. For omw use 2 iterations. All files and channels will be deconvolved with the same number of iterations.

```config.DeconIter = 20;```

Wiener filter parameter for OMW deconvolution method alpha parameter should be adjusted based on SNR and data quality. Typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7. This will take some tuning.

```config.wienerAlpha = 0.05;```

Z step size (dz). Change to step size of acquisition.

```config.dz = 0.5;```

Choose a deconvolution method. Either 'omw' or the standard matlab richardson lucy 'simplified'.

```config.RLmethod = 'simplified';```

The rest of the settings in the config struct have been set by CAMDU staff at warwick for 3i LLSM data, if you are using a different instrument it is likely you will have to change more of these settings, which can be found in the config function.

**TODO: change config to be a yaml file and provide users with a 3i or zeiss template or default then when you set common settings, select the yaml config**

**TODO: save all settings/configuration as an audit trace**

# Original PetaKit5D readme 

This project is simply a wrapper for PetaKit5D, developed by x, please cite their original work also.

# PetaKit5D

Tools for efficient and scalable processing of petabyte-scale 5D live images or large specimen images from lattice light-sheet microscopy (LLSM) and other light sheet microscopies, as well as other imaging modalities. It is featured by fast image readers and writers (for Tiff and Zarr), combined image deskew/rotation, instantly converged Richardson-Lucy (RL) deconvolution, and scalable Zarr-based stitching. It also contains some other useful tools, including 3D point detection and point tracking (based on Aguet, Upadhyayula et al., 2016), cropping, resampling, max projection, PSF analysis and visualization, and more.

## Hardware
The software works best on a Slurm-based Linux computing cluster with multiple CPU and GPU nodes for scalable, large-scale processing. For CPU nodes, we recommend at least 16 GB of RAM per core. For GPU nodes, we currently only support NVIDIA GPUs and recommend at least 24 GB of VRAM per GPU.

The software can also run on a single workstation for smaller-scale image processing tasks (up to ~1 TB, otherwise may take very long). We recommend at least 256 GB of RAM, although a smaller size may still work depending on data sizes and specific processing tasks. For relative large datasets that cause memory issues, we suggest converting the images to Zarr format and applying large-scale processing strategies.

## Usage

The tools have been tested with MATLAB R2022b-R2023a for Linux (Ubuntu 22.04), Windows (10 and 11), and MacOS (14). Toolboxes required:

`Image Processing Toolbox, Optimization Toolbox, Parallel Computing Toolbox, Signal Processing Toolbox, and Statistics and Machine Learning Toolbox.`

 Here are the steps to use the software:
1. Get the source code by either cloning the GitHub repository or downloading the ZIP file. If downloading the zip file, unzip the file to a directory.
2. Launch MATLAB, navigate to the software's root directory, and add the software to the path with `setup.m` in the command window.
````
   setup
````
3. Create a script or function to set up the workflows for your image processing tasks by calling related functions. You may follow the examples in the [demos](https://github.com/abcucberkeley/PetaKit5D/tree/main/demos). The documentation of the parameters can refer to [major_functions_documentation.txt](https://github.com/abcucberkeley/PetaKit5D/blob/main/major_functions_documentation.txt), the [GUI wiki page](https://github.com/abcucberkeley/PetaKit5D-GUI/wiki), or the parameter list in the related functions.


## Demos
The main demos for the paper:
- `demo_generic_computing_framework.m`: demo to illustrate how to use generic computing framework for user-defined functions.
- `demo_fast_tiff_zarr_readers_writers.m`: demo to illustrate Cpp-Tiff and Cpp-zarr readers and writers and compare with conventional readers and writers.
- `demo_geometric_transformation.m`: demo to illustrate how to run deskew/rotation and compare between separated and combined deskew/rotation.
- `demo_RL_deconvolution.m`: demo to illustrate how to run deconvolution and compare between traditional RL (Biggs version), WB (Guo et al. 2020), and OMW methods.
- `demo_zarr_stitching.m`: demo to illustrate how to run stitching in both skewed and DSR spaces, along with the documentation for the setup of BigStitcher (Spark version) and Stitching-Spark for the stitching benchmarks in the paper.
- `demo_large_scale_processing.m`: demo to illustrate how to set up large-scale processing for stitching, deconvolution, and deskew/rotation.
- `demo_useful_tools.m`: demo to illustrate how to set up the running for a sets of commonly used tools, i.e., resampling, cropping, max projection, tiff/zarr conversion, imaris file conversion and so on. 
- `demo_phase_and_2photon_stitching.m`: demo to illustrate how to set up the running for image list generation and stitching for 2D phase and 3D 2-photon data. 
- `demo_widefield_and_confocal_deconvolution.m`: demo to illustrate how to set up the running of deconvolution for widefield and confocal data with OMW method. 


## Python wrappers
We created Python wrappers for the main functions by calling MATLAB runtime in another repository: [PyPetaKit5D](https://github.com/abcucberkeley/PyPetaKit5D). This library is only available on Linux for now. The package of the Python wrappers can be installed using this command:
````
pip install --no-binary :all: --no-cache-dir PyPetaKit5D
````
Please refer to the [example notebooks](https://github.com/abcucberkeley/PyPetaKit5D/blob/main/notebooks) for step-by-step usage examples using our demo datasets. The parameters for the Python wrappers are identical to those in MATLAB. The documentation of the parameters can refer to [major_functions_documentation.txt](https://github.com/abcucberkeley/PetaKit5D/blob/main/major_functions_documentation.txt) or the [GUI wiki page](https://github.com/abcucberkeley/PetaKit5D-GUI/wiki).

## GUI
The software has an easy-to-use Graphical User Interface (GUI) without writing any code in another repository: [PetaKit5D-GUI](https://github.com/abcucberkeley/PetaKit5D-GUI). The GUI supports Windows, MacOS, and Linux (Ubuntu). For instructions on the installation and the usage of the PetaKit5D-GUI, visit the [GUI wiki page](https://github.com/abcucberkeley/PetaKit5D-GUI/wiki).


## Fast Tiff/Zarr readers and writers
We created independent repositories for Tiff/Zarr readers and writers for users who only need those functions: [Cpp-Tiff](https://github.com/abcucberkeley/cpp-tiff) and [Cpp-zarr](https://github.com/abcucberkeley/cpp-zarr).

We combined our Tiff/Zarr readers and [ImarisWriter](https://github.com/imaris/ImarisWriter) to develop a parallelized Imaris converter for Tiff and Zarr files. The C++ source code for this tool is available as an indepenent repository: [Parallel_Imaris_Writer](https://github.com/abcucberkeley/Parallel_Imaris_Writer).

Based on these readers and writers, we also developed a Fiji Plugin for faster reading and writing of Tiff and Zarr files within Fiji, which can be accessed from [Parallel_Fiji_Visualizer](https://github.com/abcucberkeley/Parallel_Fiji_Visualizer).


## Reference:
Please cite our paper if you find the software useful for your research:

`Ruan, X., Mueller, M., Liu, G., Görlitz, F., Fu, T., Milkie, D.E., Lillvis, J.L., Kuhn, A., Chong, J.G., Hong, J.L., Herr, C.Y.A., Hercule, W., Nienhaus, M., Killilea, A.N., Betzig, E. and Upadhyayula, S. Image processing tools for petabyte-scale light sheet microscopy data. Nature Methods (2024). https://doi.org/10.1038/s41592-024-02475-4`
