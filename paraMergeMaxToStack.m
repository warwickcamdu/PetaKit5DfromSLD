function finalStack=  paraMergeMaxToStack(inputFolder, outputFilePath, xySpacing, frameInterval, skewDirection)
    % addpath('E:\Scott\Software\Fiji.app\scripts')
    % ImageJ;
    % Ensure the output file path is a string

    if ~ischar(outputFilePath)
        outputFilePath = char(outputFilePath);
    end
    
    % List all TIFF files in the input folder
    fprintf('Listing all TIFF files in the input folder...\n');
    filePattern = fullfile(inputFolder, '*.tif');
    tiffFiles = dir(filePattern);
    
    if isempty(tiffFiles)
        error('No TIFF files found in the specified folder.');
    end
    
    % Pre-read filenames to determine the dimensions of the stack
    fprintf('Pre-reading filenames to determine stack dimensions...\n');
    timePoints = [];
    channels = [];
    
    for k = 1:length(tiffFiles)
        baseFileName = tiffFiles(k).name;
        if k == 1
            fullFileName = fullfile(inputFolder, baseFileName);
            image = readtiff_parallel(fullFileName);
            if (nargin > 4) && (skewDirection == 'Y')
                image = rot90(image,1);
            end
            imshape = size(image);
            stackSizeX = imshape(2);
            stackSizeY = imshape(1);

            
        end
        
        
        % Extract time point and channel from the filename
        tokens = regexp(baseFileName, '.*_T(\d+)_Ch(\d+).*.tif', 'tokens');
        if isempty(tokens)
            error('Filename format does not match the expected pattern.');
        end
        
        timePoint = str2double(tokens{1}{1});
        channel = str2double(tokens{1}{2});
        
        timePoints = [timePoints, timePoint];
        channels = [channels, channel];
    end
    
    % Determine the size and shape of the stack
    numTimePoints = max(timePoints) + 1
    numChannels = max(channels) + 1
    
    
    fprintf('Stack dimensions determined:\n');
    fprintf('X: %d, YY: %d, Ch: %d, Timepoints: %d\n', ...
        stackSizeX, stackSizeY, numChannels, numTimePoints);
    
    % Initialize the final stack
    fprintf('Initializing the final stack...\n');
    finalStack = zeros(stackSizeY, stackSizeX, numChannels, numTimePoints, 'uint16');

    
    fprintf('Reading images and populating the stack using parallel processing...\n');
    for k = 1:length(tiffFiles)
        baseFileName = tiffFiles(k).name
        fullFileName = fullfile(inputFolder, baseFileName)
        
        % Extract time point and channel from the filename
        tokens = regexp(baseFileName, '.*_T(\d+)_Ch(\d+).*.tif', 'tokens');
        if isempty(tokens)
            error('Filename format does not match the expected pattern.');
        end
        timePoint = str2double(tokens{1}{1});
        channel = str2double(tokens{1}{2});
        
        % Read the image
        % image = readtiff_parallel(fullFileName);
        
        % Populate the final stack
        if (nargin > 4) && (skewDirection == 'Y')
            plane = readtiff_parallel(fullFileName);
            finalStack(:, :, channel + 1, timePoint + 1) = rot90(plane,1);
        else
            finalStack(:, :, channel + 1, timePoint + 1) = readtiff_parallel(fullFileName);
        end
    end
    size(finalStack);
    
    tic;



    sizeData = [stackSizeY, stackSizeX, numChannels, numTimePoints];
    t = Tiff(outputFilePath, 'w8'); % Use 'w8' for BigTIFF
    
    % Set up the tags for the TIFF file
    tagstruct.ImageLength = sizeData(1);
    tagstruct.ImageWidth = sizeData(2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16; % Assuming 16-bit data
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter; % Use centimeters as the unit
    tagstruct.XResolution = 10000 / xySpacing; % Convert micrometers to resolution (cm)
    tagstruct.YResolution = 10000 / xySpacing; % Convert micrometers to resolution (cm)
    tagstruct.ImageDescription = sprintf(['ImageJ=1.52p\nimages=%d\nchannels=%d\nframes=%d\n' ...
                                      'hyperstack=true\nmode=grayscale\nloop=false\nfinterval=%f\n'], ...
                                      sizeData(3) * sizeData(4), sizeData(3), sizeData(4), frameInterval);

    % Write each slice
    
    for tIndex = 1:sizeData(4)
        
        
        for cIndex = 1:sizeData(3)
            % Set metadata for the current slice

            t.setTag(tagstruct);
            t.write(finalStack(:, :, cIndex, tIndex));
            t.writeDirectory();
        end
            
        
    end
    
    t.close();

    
    elapsedTime = toc;
    fprintf('5D stack saved successfully with pixel spacing and frame interval in %.2f seconds.\n', elapsedTime);
end
