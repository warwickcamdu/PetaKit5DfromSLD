% Devonvolves slidebook files (.sld) using PetaKit5D (formely known as
% llsm5dtools), followed by deskewing, and finally saving as .tif files and maximum intensity projections. 
% Petakit5d doesn't accept .sld files so firstly the sld files are converted to .tif. 

% Requires installing PetaKit5D (not the GUI version) and adding it to the matlab path 
% Installation instructions with the required matlab toolboxes are on their github:
% https://github.com/abcucberkeley/PetaKit5D

% Requires a modified Matlab bioformats toolbox to read in .sld files and
% adding it to the Matlab path.
% Possibly requires removing the bioformats that's included with petakit5d?
% This will be available somewhere like the CAMDU github

% Requires a PSF for each channel in .tif format (not .tiff).
% All the .sld files in that folder will be processed with the same PSF.

% Folder containing the .sld files to be processed.
% Don't use "C0" or "C1" in the .sld filenames or anywhere in the pathname,
% otherwise it'll break. Folder path needs to end in \


%inputFolder = 'Z:\Shared243\sbrooks\2024-06-18\to-be-deconvolvednext\';
% inputFolder = 'E:\Scott\Software\petakit5d\test-data\Series0-1_T0-1_twochannels\';
inputFolder = 'Z:\Shared243\sbrooks\2024-06-18\DeconPeta\';
% inputFolder = 'E:\Scott\Software\petakit5d\test-data\T0-2_twochannels\';
% Name of the PSF files.
% Must be .tif format and placed in the same folder as the .sld files.
% The PSF must have the same slice spacing as the image (e.g. 0.5um). 
% The metadata probably needs to be correct for the XYZ pixel spacing (e.g. 0.104 um for XY and 0.5 um for Z). 
PSF_C0 = '488_PSF.tif';
PSF_C1 = '640_PSF.tif';

% z step size
dz = 0.5;

%disable MIPS after decon, only want them after deskew
% can we output to a different directory to the tifs?
% can we delete the intermediate tifs?

% For 2024a it now uses the GPU, had to update graphics driver for matlab
% to recognise GPU, this takes a lot of pressure off of the 
% CPU, The decon is still relatively slow but has low GPU utilisation, if
% we calculate the size of each slice and how much expected GPU we can cut
% up the time points in order to parallelise the process on the GPU, expect
% 10-20x speed up so worth the time

% Choose a deconvolution method. Either 'omw' or the standard matlab richardson lucy 'simplified'. 
RLmethod = 'simplified';
% number of iterations for deconvolution. For omw use 2 iterations.
DeconIter = 20;
% Wiener filter parameter for OMW deconvolution method
% alpha parameter should be adjusted based on SNR and data quality.
% typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7
wienerAlpha = 0.05;

% Delete the raw .tif (i.e. the ones that aren't deconvolved or deskewed
deleteRawTif = false; 
% Delete the .tif files that are deconvolved but not deskewed
deleteDeconTif = false;


%% Preset Parameters 


% Deconvolution parameters 

% add the software to the path not working 
% setup([]);

% xy pixel size in um. Always 0.104 um for 3i LLSM (different to others)
xyPixelSize = 0.104;

% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.5;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'C0', 'C1', ...
                   };  

% psf path
psf_rt = inputFolder;            
PSFFullpaths = {
                [psf_rt, PSF_C0], ...
                [psf_rt, PSF_C1], ...
                };             

% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
resultDirName = 'deconvolved';

% background to subtract
Background = 100;

% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 0;
% save as 16bit; if false, save to single
Save16bit = true;
% use zarr file as input; if false, use tiff as input
zarrFile = false;
% save output as zarr file; if false,s ave as tiff
saveZarr = false;
% number of cpu cores
cpusPerTask = 4;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = false;
% use GPU for deconvolution
GPUJob = true;
% if true, save intermediate results every 5 iterations.
debug = false;
% config file for the master jobs that runs on CPU node
ConfigFile = '';
% config file for the GPU job scheduling on GPU node
GPUConfigFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;


% Deskew parameters


% also do coverslip correction rotation (usually at Warwick we don't do this)
rotate = false;
% skew angle, this is 32.8 for the 3i LLSM which is different to others
skewAngle = 32.8;
radians = deg2rad(skewAngle);

% Calculate the sine of the angle in radians
sine_value = sin(radians);
% flipZstack, I think we want this true
flipZstack = true;
% not sure this is necessary when we aren't rotating
DSRCombined = false;
% true if input is in Zarr format
zarrFile = false;
% true if saving result as Zarr files
saveZarr = false;
% true if saving result as Uint16
Save16bit = true;
% save intermediate iteration results (only for simplified, not for omw)
saveStep = false;

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;
% use master job for task computing or not. 
masterCompute = true;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;


%% Step 1. Convert the .sld files into .tif files

% Find the .sld files in the directory
filePattern = fullfile(inputFolder, '*.sld'); 
theFiles = dir(filePattern);

% Store total number of sld files to study
nFiles = length(theFiles);

% Iterate through all the .sld files in the directory
for k = 1:nFiles
    fprintf("   >> Converting .sld to tif: %3d / %3d\n", k, nFiles);

    % Define full file name for current loop iteration
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName)

    r = bfGetReader(fullFileName);

    %access the OME metadata and get number of series
    omeMeta = r.getMetadataStore();
    nSeries = r.getSeriesCount();

    %Iterate through series within the file
    for S = 0:nSeries-1


        %switch between series and load that series
        r.setSeries(S);
        %r.getSeries();

        %get metadata and extract important features
        % check X and Y are correct and not switched
        omeMeta = r.getMetadataStore();
        stackSizeX = omeMeta.getPixelsSizeX(0).getValue();      %image width in pixels
        stackSizeY = omeMeta.getPixelsSizeY(0).getValue();     %image height in pixels
        stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue();      %number of slices
        stackSizeC = omeMeta.getPixelsSizeC(0).getValue();      %number of channels
        stackSizeT = omeMeta.getPixelsSizeT(0).getValue();      %number of time points
        % Extract physical pixel size (XY spacing)
        pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0); % in micrometers
        if ~isempty(pixelSizeX)
            pixelSizeX = double(pixelSizeX.value());  
        else
            pixelSizeX = NaN;
        end

        pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0); % in micrometers
        if ~isempty(pixelSizeY)
            pixelSizeY = double(pixelSizeY.value());
        else
            pixelSizeY = NaN;
        end

        % Extract Z spacing
        pixelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0); % in micrometers
        if ~isempty(pixelSizeZ)
            pixelSizeZ = double(pixelSizeZ.value());
        else
            pixelSizeZ = NaN;
        end
        deskewedZSpacing = sine_value * pixelSizeZ;

        % Extract frame interval
        frameInterval = 0;
        % % Extract frame interval
        % frameInterval = omeMeta.getPixelsTimeIncrement(1); % in seconds
        % if ~isempty(frameInterval)
        %     frameInterval = frameInterval.value();
        % else
        %     frameInterval = NaN;
        % end
        if stackSizeC == 1
            PSFFullpaths = {
                [psf_rt, PSF_C0], ...
                };     
            ChannelPatterns = {'C0', ...
                   };  
        end 



        % Print extracted values
        fprintf('Stack Size (X, Y, Z, C, T): (%d, %d, %d, %d, %d)\n', stackSizeX, stackSizeY, stackSizeZ, stackSizeC, stackSizeT);
        fprintf('Pixel Size (X, Y): (%.3f, %.3f) micrometers\n', pixelSizeX, pixelSizeY);
        fprintf('Z Spacing: %.2f micrometers\n', pixelSizeZ);
        fprintf('Deskewed Z Spacing: %.3f micrometers\n', deskewedZSpacing);


        seriesName = char(omeMeta.getImageName(S));

        % Take the .sld filename and the series image name and make a new
        % folder based on this
        seriesFolderName = strrep(baseFileName, ".sld", "");
        seriesNameNoSpaces = strrep(seriesName, " ", "_");
        currentSeriesFolder = seriesFolderName+'_'+seriesNameNoSpaces;
        mkdir(inputFolder,currentSeriesFolder);

        currentSeriesPath = fullfile(inputFolder, currentSeriesFolder);

        % make a folder to store the .tif files


        mkdir(currentSeriesPath,'tifs');
        tifDir = currentSeriesPath+ '\'+ 'tifs';

        %skip the series if it has only one Z-slice
        plane_count = 0;
        if stackSizeZ>1
            %We need to store all of Z-stacks of this time-point and
            %channel in an array to be processed later, so set up and empty array
            %and start a count
            count = 1;
            array = [];

            %iterate through all the timepoints
            for T = 0:stackSizeT-1
                %iterate through all the channels
                for C = 0:stackSizeC-1
                    %iterate through all the z-slices
                    for Z = 0:stackSizeZ-1

                        %Use the index to read in the specific plane and
                        %convert to double
                        plane = bfGetPlane(r, r.getIndex(Z, C, T) +1);
                        plane = double(plane);

                        %Add plane to array at position (count, 1)(in essence
                        %you are appending the array) and add 1 to count.
                        %Possibly this step could be improved, there's
                        %maybe a simpler way to read this without a for
                        %loop
                        array(:,:,count) = plane;
                        count = count+1;
                        plane_count = plane_count+1;

                        % this should be first of second timepoint
                        if plane_count == int32(stackSizeZ*stackSizeC)

                            frameInterval = omeMeta.getPlaneDeltaT(S, plane_count).value().doubleValue()/1000; % in seconds
                            firstframeInterval = omeMeta.getPlaneDeltaT(0, 0).value().doubleValue()/1000; % in seconds
                            frameInterval =  frameInterval - firstframeInterval;
                        end

                    end

                    %prepare the stack array to be saved
                    outputArray = array(1:stackSizeY, 1:stackSizeX, 1:stackSizeZ);
                    outputArray = uint16(outputArray);

                    %prepare the image name
                    strSld = baseFileName(1:end-4);
                    strS = num2str(S);
                    strT = num2str(T); 
                    strT = pad(strT,4,'left','0'); % zero padding
                    strC = num2str(C);

                    % Ensure all variables are character arrays
                    tifDir = char(tifDir);
                    tifFullpath = [tifDir '\' strSld '_S' strS '_T' strT '_C' strC '.tif'];

                    %save the array as a tif
                    %I think this doesn't save metadata but that doesn't
                    %seem to matter for the deconvolution step
                    parallelWriteTiff(tifFullpath,outputArray);

                    %clear array for next channel
                    array = [];
                    count = 1;
                end

                %clear array for next timepoint
                array = [];
                count = 1;


            end
        end
        fprintf('Frame Interval: %.2f seconds\n', frameInterval);

        % After we have saved all the individual 3D stacks for each image
        % the current metadata will be correct and we can run deconv and
        % deskew as we go along then reconstruct into a new 5D image

        %% Step 2: Deconvolution 
        %% Step 2.1: set parameters 
        % set in top section
        
        %% Step 2.2: run the deconvolution with given parameters. 
        % the results will be saved in matlab_decon under the dataPaths. 
        % the next step is deskew/rotate (if in skewed space for x-stage scan) or 
        % rotate (if objective scan) or other processings. 
        fprintf('Starting deconvolution...\n\n');
        
        
        XR_decon_data_wrapper(tifDir, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
            'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
            'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
            'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
            'Background', Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', DeconIter, ...
            'fixIter', fixIter, 'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, ...
            'zarrFile', zarrFile, 'saveZarr', saveZarr, 'parseCluster', parseCluster, ...
            'largeFile', largeFile, 'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask, ...
            'ConfigFile', ConfigFile, 'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);
        
        % release GPU if using GPU computing
        if GPUJob && gpuDeviceCount('available') > 0
            reset(gpuDevice);
        end
        
        %% Step 3: deskew the deconvolved results
        
        dataPath_exps = [tifDir '\' resultDirName];
        
        
        
        XR_deskew_rotate_data_wrapper(dataPath_exps, skewAngle=skewAngle, flipZstack=flipZstack, DSRCombined=DSRCombined, rotate=rotate, xyPixelSize=xyPixelSize, dz=dz, ...
            Reverse=Reverse, ChannelPatterns=ChannelPatterns, largeFile=largeFile, ...
            zarrFile=zarrFile, saveZarr=saveZarr, Save16bit=Save16bit, parseCluster=parseCluster, ...
            masterCompute=masterCompute, configFile=configFile, mccMode=mccMode);
        
        %outputTiffFile = currentSeriesFolder + ".tif";
        outputTiffFile = currentSeriesPath + ".tif";
        outputTiffPath = fullfile(inputFolder, outputTiffFile);
        inputToMerge = [dataPath_exps '\' 'DS'];
        paraMergeTiffFilesToMultiDimStack(inputToMerge, outputTiffFile,pixelSizeX, deskewedZSpacing, frameInterval);
        
        
    end

end




% This was not clear, may need to be added just before the deconvolution
% step
% % move to the PetaKit5D root directory
% curPath = pwd;
% if ~endsWith(curPath, 'PetaKit5D')
%     mfilePath = mfilename('fullpath');
%     if contains(mfilePath,'LiveEditorEvaluationHelper')
%         mfilePath = matlab.desktop.editor.getActiveFilename;
%     end
% 
%     mPath = fileparts(mfilePath);
%     if endsWith(mPath, 'demos')
%         cd(mPath);
%         cd('..')
%     end
% end



%% Step 4: delete intermediate .tif files
 
% deletes the raw tifs if the flag is true
if deleteRawTif == true 

    % Find the raw .tif files (i.e. not deconvolved, not deskewed)
    filePattern = fullfile(tifDir, '*.tif'); 
    theFiles = dir(filePattern);
    
    % Store total number of sld files to study
    nFiles = length(theFiles);
    
    % Iterate through all the .sld files in the directory
    for k = 1:nFiles        
    
        % Define full file name for current loop iteration
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        delete(fullFileName);
    end
end


% deletes the .tif files that are deconvolved but not deskewed, if the flag is true
if deleteDeconTif == true 

    % Path to the deconvolved .tif files
    deconTifDir = [tifDir '\' resultDirName];
    % Find the raw .tif files (i.e. not deconvolved, not deskewed)
    
    filePattern = fullfile(deconTifDir, '*.tif'); 
    theFiles = dir(filePattern);
    
    % Store total number of sld files to study
    nFiles = length(theFiles);
    
    % Iterate through all the .sld files in the directory
    for k = 1:nFiles
            
        % Define full file name for current loop iteration
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        delete(fullFileName);
    end
end


