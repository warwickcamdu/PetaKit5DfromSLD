function finalStack=  paraMergeTiffFilesToMultiDimStack(inputFolder, outputFilePath, xySpacing, zSpacing, frameInterval)
    % addpath('E:\Scott\Software\Fiji.app\scripts')
    % ImageJ;
    % Ensure the output file path is a string
    % TODO test with 1 channel 1 timepoint
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
            imshape = size(image);
            stackSizeX = imshape(2);
            stackSizeY = imshape(1);
            stackSizeZ = imshape(3);

            
        end
        
        
        % Extract time point and channel from the filename
        tokens = regexp(baseFileName, '.*_T(\d+)_Ch(\d+).tif', 'tokens');
        if isempty(tokens)
            error('Filename format does not match the expected pattern.');
        end
        
        timePoint = str2double(tokens{1}{1});
        channel = str2double(tokens{1}{2});
        
        timePoints = [timePoints, timePoint];
        channels = [channels, channel];
    end
    
    % Determine the size and shape of the stack
    numTimePoints = max(timePoints) + 1;
    numChannels = max(channels) + 1;
    
    
    fprintf('Stack dimensions determined:\n');
    fprintf('X: %d, YY: %d, Z: %d, Ch: %d, Timepoints: %d\n', ...
        stackSizeX, stackSizeY, stackSizeZ, numChannels, numTimePoints);
    
    % Initialize the final stack
    fprintf('Initializing the final stack...\n');
    finalStack = zeros(stackSizeY, stackSizeX, stackSizeZ, numChannels, numTimePoints, 'uint16');
    numTimePoints;

    
    fprintf('Reading images and populating the stack using parallel processing...\n');
    for k = 1:length(tiffFiles)
        baseFileName = tiffFiles(k).name;
        fullFileName = fullfile(inputFolder, baseFileName);
        
        % Extract time point and channel from the filename
        tokens = regexp(baseFileName, '.*_T(\d+)_Ch(\d+).tif', 'tokens');
        if isempty(tokens)
            error('Filename format does not match the expected pattern.');
        end
        timePoint = str2double(tokens{1}{1});
        channel = str2double(tokens{1}{2});
        
        % Read the image
        % image = readtiff_parallel(fullFileName);
        
        % Populate the final stack
        finalStack(:, :, :, channel + 1, timePoint + 1) = readtiff_parallel(fullFileName);
    end
    % Extract the slice
    % slice = finalStack(:, :, 45, 2, 1);
    % 
    % % Normalize the slice
    % sliceNormalized = mat2gray(slice);
    % 
    % % Display the normalized slice
    % imshow(sliceNormalized);
    % title('Normalized Slice');

    
    
    % Verify metadata
    % zSpacing = 0.271; 
    % frameInterval = 4;
    % xySpacing = 0.104;

    
    
    
    % % Save the final stack as an OME-TIFF using Bio-Formats
    % fprintf('Saving the final stack...\n');
    % tic;
    % bfsave(finalStack, outputFilePath);
    % elapsedTime = toc;
    % fprintf('5D stack saved successfully with pixel spacing and frame interval in %.2f seconds.\n', elapsedTime);
    % IJ.saveAsTiff(finalStack, outputFilePath);

    % Assuming yourData is your 5D stack in YXZCT order
    tic;

    sizeData = [stackSizeY, stackSizeX, stackSizeZ, numChannels, numTimePoints];
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
    tagstruct.ImageDescription = sprintf(['ImageJ=1.52p\nimages=%d\nchannels=%d\nslices=%d\nframes=%d\n' ...
                                      'hyperstack=true\nmode=grayscale\nloop=false\nspacing=%f\n' ...
                                      'unit=micron\nfinterval=%f\n'], ...
                                      sizeData(3) * sizeData(4) * sizeData(5), sizeData(4), sizeData(3), sizeData(5), zSpacing, frameInterval);

    % Write each slice
    
    for tIndex = 1:sizeData(5)
        
        for zIndex = 1:sizeData(3)
            for cIndex = 1:sizeData(4)
                % Set metadata for the current slice

                t.setTag(tagstruct);
                t.write(finalStack(:, :, zIndex, cIndex, tIndex));
                t.writeDirectory();
            end
            
        end
    end
    
    t.close();

    
    elapsedTime = toc;
    fprintf('5D stack saved successfully with pixel spacing and frame interval in %.2f seconds.\n', elapsedTime);
end
