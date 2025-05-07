function modularPipeline()
% modularPipeline prompts the user for:
%    1. A folder containing PSF files. Assumes PSF filenames contain the text
%       'PSF_CHX' (where X is a number) and uses those to build the PSF paths.
%    2. A folder containing the files to deconvolve (either .sld or .tif).
%
% The code then processes the data using one or both pipelines (decon+deskew and/or
% deskew-only), with z–axis padding applied during conversion and removed before merging.

    %% --- UI: Ask the User for Required Folders ---
    psfFolder = uigetdir([], 'Select the folder containing your PSF files');
    if psfFolder == 0
        error('No PSF folder selected.');
    end
    
    inputFolder = uigetdir([], 'Select the folder containing files to deconvolve (SLD, CZI or TIFF)');
    if inputFolder == 0
        error('No input folder selected.');
    end
    
    %% --- Get Default Configuration and Update with User Choices ---
    config = getDefaultConfig();
    config.inputFolder = inputFolder;
    
    % Build the PSF file list by scanning the PSF folder.
    psfFiles = dir(fullfile(psfFolder, '*PSF_CH*.tif'));
    if isempty(psfFiles)
        error('No PSF files found in %s', psfFolder);
    end
    
    % Extract the channel number from each filename.
    psfArray = struct('channel', {}, 'fullpath', {});
    for i = 1:length(psfFiles)
        fname = psfFiles(i).name;
        tokens = regexp(fname, 'PSF_CH(\d+)', 'tokens');
        if ~isempty(tokens)
            channelNum = str2double(tokens{1}{1});
            psfArray(end+1).channel = channelNum;  %#ok<AGROW>
            psfArray(end).fullpath = fullfile(psfFiles(i).folder, fname);
        end
    end
    
    if isempty(psfArray)
        error('No valid PSF files found with pattern "PSF_CHX" in %s', psfFolder);
    end
    
    % Sort the PSF files by channel number.
    [~, idx] = sort([psfArray.channel]);
    psfArray = psfArray(idx);
    
    % Update the configuration with PSF file paths and channel patterns.
    config.PSFFullpaths = cell(1, length(psfArray));
    config.ChannelPatterns = cell(1, length(psfArray));
    for i = 1:length(psfArray)
        config.PSFFullpaths{i} = psfArray(i).fullpath;
        % Create a channel pattern (e.g. 'Ch1', 'Ch2', ...).
        config.ChannelPatterns{i} = ['Ch' num2str(psfArray(i).channel)];
    end

    fprintf('Selected PSF files:\n');
    disp(config.PSFFullpaths);
    fprintf('Channel Patterns:\n');
    disp(config.ChannelPatterns);
    fprintf('Processing files in folder: %s\n', config.inputFolder);
    
    %% --- Process the Input Data ---
    % Determine if the input folder contains .sld, .czi or .tif files.
    sldyFiles = dir(fullfile(config.inputFolder, '*.sldy'));
    sldFiles = dir(fullfile(config.inputFolder, '*.sld'));
    allTifFiles = dir(fullfile(config.inputFolder, '*.tif'));
    allTifFiles = allTifFiles(~[allTifFiles.isdir]);
    cziFiles = dir(fullfile(config.inputFolder, '*.czi'));
    
    if ~isempty(sldFiles)
         % Process each SLD file.
         for i = 1:length(sldFiles)
             sldFullPath = fullfile(sldFiles(i).folder, sldFiles(i).name);
             processSldFile(sldFullPath, config);
         end
    elseif ~isempty(cziFiles)
        %default czi config
        config = getCziDefaultConfig(config);
        %Process each CZI file
        for i=1:length(cziFiles)
            cziFullPath = fullfile(cziFiles(i).folder, cziFiles(i).name);
            processSldFile(cziFullPath, config);
        end
    elseif ~isempty(sldyFiles)
        %Process each SLDY file
        for i=1:length(sldyFiles)
            sldyFullPath = fullfile(sldyFiles(i).folder, sldyFiles(i).name);
            processSldFile(sldyFullPath, config);
        end
    elseif ~isempty(allTifFiles)
         % Pattern to match _T<number>_Ch<number>
         pattern = '_T\d+_Ch\d+';

         % Filter TIF files that match the pattern
         tifFiles = [];
         for i = 1:length(allTifFiles)
            filename = allTifFiles(i).name;
            if ~isempty(regexp(filename, pattern, 'once'))
                tifFiles = [tifFiles; allTifFiles(i)];  % Append matching file
            end
         end

        if ~isempty(tifFiles)
            % Only process matching TIF files
            filePaths = fullfile({tifFiles.folder}, {tifFiles.name});
            disp('Processing series of 3D TIF files...');
            seriesResult = processTifFolder(config);
            if isempty(seriesResult)
                 error('No valid TIFF series found in %s', config.inputFolder);
            else
                if strcmp(config.processingMode, 'decon+deskew') || strcmp(config.processingMode, 'both')
                    runDeconDeskewPipeline(seriesResult, config);
                end
                if strcmp(config.processingMode, 'deskew-only') || strcmp(config.processingMode, 'both')
                    runDeskewOnlyPipeline(seriesResult, config);
                end
                deleteIntermediateFiles(seriesResult.tifDir, config);
            end
        else
            % No matching files, do something else
            disp('Processing Tif files...');
            for i = 1:length(allTifFiles)
                filepath = fullfile(allTifFiles.folder,allTifFiles(i).name);
                processSldFile(filepath, config);
            end
        end
    else
         error('No .sld or .tif files found in folder %s', config.inputFolder);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: processSldFile
function processSldFile(sldFileName, config)
    % Open the .sld file using Bio-Formats.
    r = bfGetReader(sldFileName);
    omeMeta = r.getMetadataStore();
    if endsWith(sldFileName, {'.sld','.sldy'}, 'IgnoreCase', true)
        nSeries = r.getSeriesCount();
    else
        nSeries=1;
    end
    
    % Process each series in the .sld file.
    for S = 0:nSeries-1
        if endsWith(sldFileName, {'.sld','.sldy'}, 'IgnoreCase', true)
        r.setSeries(S);
        end
        seriesResult = convertSeriesToTif(r, S, sldFileName, config);
        if isempty(seriesResult)
            continue;  % Skip series with only one Z-slice.
        end
        
        if strcmp(config.processingMode, 'decon+deskew') || strcmp(config.processingMode, 'both')
            runDeconDeskewPipeline(seriesResult, config);
        end
        if strcmp(config.processingMode, 'deskew-only') || strcmp(config.processingMode, 'both')
            runDeskewOnlyPipeline(seriesResult, config);
        end
        deleteIntermediateFiles(seriesResult.tifDir, config);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: convertSeriesToTif
function seriesResult = convertSeriesToTif(r, seriesIndex, sldFileName, config)
    omeMeta = r.getMetadataStore();
    % Extract image dimensions.
    stackSizeX = omeMeta.getPixelsSizeX(seriesIndex).getValue();
    stackSizeY = omeMeta.getPixelsSizeY(seriesIndex).getValue();
    stackSizeZ = omeMeta.getPixelsSizeZ(seriesIndex).getValue();
    stackSizeC = omeMeta.getPixelsSizeC(seriesIndex).getValue();
    stackSizeT = omeMeta.getPixelsSizeT(seriesIndex).getValue();
    
    % Get physical pixel sizes. If unavailable (or NaN) then substitute defaults from config.
    pixelSizeX_obj = omeMeta.getPixelsPhysicalSizeX(seriesIndex);
    if isempty(pixelSizeX_obj)
        pixelSizeX = config.xyPixelSize;
    else
        pixelSizeX = double(pixelSizeX_obj.value());
        if isnan(pixelSizeX)
            pixelSizeX = config.xyPixelSize;
        end
    end
    
    pixelSizeY_obj = omeMeta.getPixelsPhysicalSizeY(seriesIndex);
    if isempty(pixelSizeY_obj)
        pixelSizeY = config.xyPixelSize;
    else
        pixelSizeY = double(pixelSizeY_obj.value());
        if isnan(pixelSizeY)
            pixelSizeY = config.xyPixelSize;
        end
    end
    
    pixelSizeZ_obj = omeMeta.getPixelsPhysicalSizeZ(seriesIndex);
    if isempty(pixelSizeZ_obj)
        pixelSizeZ = config.dz;
    else
        pixelSizeZ = double(pixelSizeZ_obj.value());
        if isnan(pixelSizeZ)
            pixelSizeZ = config.dz;
        end
    end
    
    % Calculate the deskewed Z spacing using the Z spacing value.
    deskewedZSpacing = sin(deg2rad(config.skewAngle)) * pixelSizeZ;
    
    frameInterval = 0;
    plane_count = 0;
    
    % Create an output folder based on the series.
    seriesName = char(omeMeta.getImageName(seriesIndex));
    [~, baseFileName, ~] = fileparts(sldFileName);
    seriesNameNoSpaces = strrep(seriesName, ' ', '_');
    currentSeriesFolder = [baseFileName, '_', seriesNameNoSpaces];
    currentSeriesPath = fullfile(config.inputFolder, currentSeriesFolder);
    if ~exist(currentSeriesPath, 'dir')
        mkdir(currentSeriesPath);
    end
    tifDir = fullfile(currentSeriesPath, 'tifs');
    if ~exist(tifDir, 'dir')
        mkdir(tifDir);
    end
    
    fprintf('Series %d: Dimensions (X,Y,Z,C,T) = (%d,%d,%d,%d,%d)\n', ...
             seriesIndex, stackSizeX, stackSizeY, stackSizeZ, stackSizeC, stackSizeT);
    fprintf('Pixel Size (X,Y): (%.3f, %.3f) um, Z Spacing: %.2f um, Deskewed Z Spacing: %.3f um\n', ...
             pixelSizeX, pixelSizeY, pixelSizeZ, deskewedZSpacing);
    
    if stackSizeZ <= 1
        warning('Skipping series %d as it contains only one Z-slice.', seriesIndex);
        seriesResult = [];
        return;
    end
    
    if endsWith(sldFileName, {'.tif', '.tiff'}, 'IgnoreCase', true)
        fullarray=readtiff_parallel(sldFileName);
    end
    % Loop through timepoints and channels.
    for T = 0:stackSizeT-1
        for C = 0:stackSizeC-1
            array = [];
            count = 1;
            
            if endsWith(sldFileName, {'.sld','.sldy','.czi'}, 'IgnoreCase', true)
                for Z = 0:stackSizeZ-1
                    plane = bfGetPlane(r, r.getIndex(Z, C, T) + 1);
                    array(:, :, count) = double(plane);
                    count = count + 1;
                    plane_count = plane_count + 1;
                    if plane_count == int32(stackSizeZ * stackSizeC) + 1
                        frameInterval = omeMeta.getPlaneDeltaT(seriesIndex, plane_count).value().doubleValue()/1000;
                        firstframeInterval = omeMeta.getPlaneDeltaT(seriesIndex, 0).value().doubleValue()/1000;
                        frameInterval = frameInterval - firstframeInterval;
                    end
                end
            elseif endsWith(sldFileName, {'.tif', '.tiff'}, 'IgnoreCase', true)
                array=fullarray(:,:,(C+1):stackSizeC:2*stackSizeZ*(T+1));
            end
            
            outputArray = uint16(array(:, :, 1:stackSizeZ));
            outputArray = applyZPadding(outputArray, config);
            
            strS = num2str(seriesIndex);
            strT = pad(num2str(T), 4, 'left', '0');
            strC = num2str(C);
            tifFileName = fullfile(tifDir, sprintf('%s_S%s_T%s_Ch%s.tif', baseFileName, strS, strT, strC));
            parallelWriteTiff(tifFileName, outputArray);
        end
    end
    
    seriesResult.tifDir = tifDir;
    seriesResult.currentSeriesFolder = currentSeriesFolder;
    seriesResult.currentSeriesPath = currentSeriesPath;
    seriesResult.frameInterval = frameInterval;
    seriesResult.pixelSizeX = pixelSizeX;
    seriesResult.pixelSizeY = pixelSizeY;
    seriesResult.deskewedZSpacing = deskewedZSpacing;
end

%% -----------------------------------------------------------------------
%% Local Function: processTifFolder
function seriesResult = processTifFolder(config)
    % For a folder of TIFF files, assume that the folder itself contains the raw images.
    tifDir = config.inputFolder;
    [~, currentSeriesFolder, ~] = fileparts(tifDir);
    if isempty(currentSeriesFolder)
       currentSeriesFolder = 'raw_tif_series';
    end
    
    seriesResult.tifDir = tifDir;
    seriesResult.currentSeriesFolder = currentSeriesFolder;
    seriesResult.currentSeriesPath = tifDir;
    seriesResult.frameInterval = 0;  % Default if metadata is unavailable.
    seriesResult.pixelSizeX = config.xyPixelSize;
    seriesResult.deskewedZSpacing = sin(deg2rad(config.skewAngle)) * config.dz;
end

%% -----------------------------------------------------------------------
%% Local Function: applyZPadding
function paddedArray = applyZPadding(array, config)
    switch config.z_edge_padding
        case 'none'
            paddedArray = array;
        case 'zero'
            paddedArray = padarray(array, [0, 0, config.z_padding], 0, 'both');
        case 'mirror'
            paddedArray = padarray(array, [0, 0, config.z_padding], 'symmetric', 'both');
        case 'gaussian'
            frontPad = config.gaussian_mean + config.gaussian_std .* randn(size(array,1), size(array,2), config.z_padding);
            backPad  = config.gaussian_mean + config.gaussian_std .* randn(size(array,1), size(array,2), config.z_padding);
            paddedArray = cat(3, frontPad, array, backPad);
        case 'fixed'
            frontPad = config.fixed_value * ones(size(array,1), size(array,2), config.z_padding);
            backPad  = config.fixed_value * ones(size(array,1), size(array,2), config.z_padding);
            paddedArray = cat(3, frontPad, array, backPad);
        otherwise
            error('Invalid z_edge_padding option: %s', config.z_edge_padding);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: runDeconDeskewPipeline
function runDeconDeskewPipeline(seriesResult, config)
    fprintf('Running deconvolution+deskew pipeline for series: %s\n', seriesResult.currentSeriesFolder);
    
    % Deconvolution step.
    XR_decon_data_wrapper(seriesResult.tifDir, 'resultDirName', config.resultDirName, 'xyPixelSize', config.xyPixelSize, ...
                'dz', config.dz, 'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'PSFFullpaths', config.PSFFullpaths, ...
                'dzPSF', config.dzPSF, 'parseSettingFile', config.parseSettingFile, 'RLmethod', config.RLmethod, ...
                'wienerAlpha', config.wienerAlpha, 'OTFCumThresh', config.OTFCumThresh, 'skewed', config.skewed, ...
                'Background', config.Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', config.DeconIter, ...
                'fixIter', config.fixIter, 'EdgeErosion', config.EdgeErosion, 'Save16bit', config.Save16bit, ...
                'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'parseCluster', config.parseCluster, ...
                'largeFile', config.largeFile, 'GPUJob', config.GPUJob, 'debug', config.debug, 'cpusPerTask', config.cpusPerTask, ...
                'ConfigFile', config.ConfigFile, 'GPUConfigFile', config.GPUConfigFile, 'mccMode', config.mccMode);
            
    if config.GPUJob && gpuDeviceCount('available') > 0
         reset(gpuDevice);
    end
    
    % Deskew step.
    dataPath_exps = fullfile(seriesResult.tifDir, config.resultDirName);
    XR_deskew_rotate_data_wrapper(dataPath_exps, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Remove z-padding from the decon+deskew results.
    deconDSDir = fullfile(dataPath_exps, 'DS');
    removePaddingFromDir(deconDSDir, config);
    
    % Merge the deconvolved+deskewed images.
    outputTiffFile = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '.tif']);
    paraMergeTiffFilesToMultiDimStack(deconDSDir, outputTiffFile, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    
    outputTiffFileMax = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_MAX.tif']);
    inputToMergeMax = fullfile(deconDSDir, 'MIPs');
    paraMergeMaxToStack(inputToMergeMax, outputTiffFileMax, seriesResult.pixelSizeX, seriesResult.frameInterval);
end

%% -----------------------------------------------------------------------
%% Local Function: runDeskewOnlyPipeline
function runDeskewOnlyPipeline(seriesResult, config)
    fprintf('Running deskew-only pipeline for series: %s\n', seriesResult.currentSeriesFolder);
    
    XR_deskew_rotate_data_wrapper(seriesResult.tifDir, 'resultDirName', config.resultDirNameDeskew, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Remove padding from the deskew-only output.
    deskewDSDir = fullfile(seriesResult.tifDir, config.resultDirNameDeskew);
    deskewDSDir
    removePaddingFromDir(deskewDSDir, config);
    
    % Merge the deskew-only images.
    outputTiffFileDeskew = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_deskew.tif']);
    paraMergeTiffFilesToMultiDimStack(deskewDSDir, outputTiffFileDeskew, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    
    outputTiffFileDeskewMax = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_deskew_MAX.tif']);
    inputToMergeDeskewMax = fullfile(deskewDSDir, 'MIPs');
    paraMergeMaxToStack(inputToMergeDeskewMax, outputTiffFileDeskewMax, seriesResult.pixelSizeX, seriesResult.frameInterval);
end

%% -----------------------------------------------------------------------
%% Local Function: removePaddingFromDir
function removePaddingFromDir(targetDir, config)
    fileList = dir(fullfile(targetDir, '*.tif'));
    if ~strcmp(config.z_edge_padding, 'none')
        for i = 1:length(fileList)
            filePath = fullfile(targetDir, fileList(i).name);
            img = parallelReadTiff(filePath);
            if size(img, 3) > 2 * config.z_padding
                img_no_padding = img(:, :, (config.z_padding+1):(end-config.z_padding));
                parallelWriteTiff(filePath, img_no_padding);
            else
                warning('Not enough depth for padding removal in file: %s', fileList(i).name);
            end
        end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: deleteIntermediateFiles
function deleteIntermediateFiles(tifDir, config)
    if config.deleteRawTif
         rawFiles = dir(fullfile(tifDir, '*.tif'));
         for k = 1:length(rawFiles)
             delete(fullfile(rawFiles(k).folder, rawFiles(k).name));
         end
    end
    if config.deleteDeconTif
         deconTifDir = fullfile(tifDir, config.resultDirName);
         deconFiles = dir(fullfile(deconTifDir, '*.tif'));
         for k = 1:length(deconFiles)
             delete(fullfile(deconFiles(k).folder, deconFiles(k).name));
         end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: getDefaultConfig
function config = getDefaultConfig()
    % Folder settings (these will be updated by the user selections).
    config.inputFolder = '';   % The folder to deconvolve (set via UI)
    
    % PSF related settings will be updated via UI:
    config.PSFFullpaths = {};  % Cell array of PSF file paths.
    config.ChannelPatterns = {}; % Patterns (e.g., 'Ch1','Ch2', ...) for each channel.
    
    % Imaging parameters.
    config.dz = 0.5;
    config.xyPixelSize = 0.104;
    
    % z–axis padding settings.
    config.z_edge_padding = 'gaussian';   % Options: 'none', 'zero', 'mirror', 'gaussian', 'fixed'
    config.z_padding = 20;
    config.gaussian_mean = 102.27;
    config.gaussian_std = 3.17;
    config.fixed_value = 100;
    
    % Flags for deleting intermediate files.
    config.deleteRawTif = false;
    config.deleteDeconTif = false;
    
    % Deconvolution parameters.
    config.RLmethod = 'simplified';
    config.DeconIter = 10;
    config.wienerAlpha = 0.05;
    
    % Acquisition and PSF parameters.
    config.Reverse = true;
    config.dzPSF = 0.5;
    config.parseSettingFile = false;
    % The ChannelPatterns and PSFFullpaths will be updated from the PSF folder.
    config.OTFCumThresh = 0.9;
    config.skewed = true;
    config.Background = 100;
    config.fixIter = true;
    config.EdgeErosion = 0;
    config.Save16bit = true;
    config.zarrFile = false;
    config.saveZarr = false;
    config.cpusPerTask = 4;
    config.parseCluster = false;
    config.largeFile = false;
    config.GPUJob = true;
    config.debug = false;
    config.ConfigFile = '';
    config.GPUConfigFile = '';
    config.mccMode = false;
    
    % Deskew and rotation parameters.
    config.rotate = false;
    config.skewAngle = 32.8;
    config.flipZstack = true;
    config.DSRCombined = false;
    config.masterCompute = true;
    config.configFile = '';
    
    % Processing mode: Choose 'deskew-only', 'decon+deskew', or 'both'.
    config.processingMode = 'both';
    
    % Output folder names.
    config.resultDirName = 'deconvolved';       % For deconvolution+deskew branch.
    config.resultDirNameDeskew = 'DS';   % For deskew-only branch.
end

function config = getCziDefaultConfig(config)
    config.xyPixelSize = 0.1449922;
    config.skewAngle = 32.8;
    config.flipZstack = false;
end