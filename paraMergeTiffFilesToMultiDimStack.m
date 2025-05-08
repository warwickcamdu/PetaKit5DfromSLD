function paraMergeTiffFilesToMultiDimStack(inputFolder, outputFilePath, xySpacing, zSpacing, frameInterval)
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
    
    % Determine the stack dimensions based on the first file
    baseFileName = tiffFiles(1).name;
    fullFileName = fullfile(inputFolder, baseFileName);
    image = readtiff_parallel(fullFileName);
    [stackSizeY, stackSizeX, stackSizeZ] = size(image);

fileMeta = struct('name', {}, 'T', {}, 'Ch', {});
for k = 1:length(tiffFiles)
    tokens = regexp(tiffFiles(k).name, '.*_T(\d+)_Ch(\d+).tif', 'tokens');
    if isempty(tokens)
        error('Filename format does not match the expected pattern: %s', tiffFiles(k).name);
    end

    timePoint = str2double(tokens{1}{1});
    channel = str2double(tokens{1}{2});

    fileMeta(end+1).name = tiffFiles(k).name; %#ok<AGROW>
    fileMeta(end).T = timePoint;
    fileMeta(end).Ch = channel;
end

    %Sort files
    [~, sortIdx] = sortrows([[fileMeta.T]', [fileMeta.Ch]']);
    fileMeta = fileMeta(sortIdx);
    % Determine the number of time points and channels
    allT = [fileMeta.T];
allCh = [fileMeta.Ch];
numTimePoints = max(allT) + 1;
numChannels = max(allCh) + 1;

    fprintf('Stack dimensions determined:\n');
    fprintf('X: %d, Y: %d, Z: %d, Ch: %d, Timepoints: %d\n', ...
        stackSizeX, stackSizeY, stackSizeZ, numChannels, numTimePoints);

    % Set up TIFF for writing
    tic;
    t = Tiff(outputFilePath, 'w8'); % Use 'w8' for BigTIFF

    % Set up the tags for the TIFF file
    tagstruct.ImageLength = stackSizeY;
    tagstruct.ImageWidth = stackSizeX;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16; % Assuming 16-bit data
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
    tagstruct.XResolution = 10000 / xySpacing;
    tagstruct.YResolution = 10000 / xySpacing;
    tagstruct.ImageDescription = sprintf(['ImageJ=1.52p\nimages=%d\nchannels=%d\nslices=%d\nframes=%d\n' ...
                                          'hyperstack=true\nmode=grayscale\nloop=false\nspacing=%f\n' ...
                                          'unit=micron\nfinterval=%f\n'], ...
                                          stackSizeZ * numChannels * numTimePoints, numChannels, stackSizeZ, numTimePoints, zSpacing, frameInterval);

    % Process each file, read and write directly to the TIFF
    fprintf('Reading images and writing to the stack plane-by-plane...\n');
    for tIndex = 1:numTimePoints
        intermediateStack = zeros(stackSizeY, stackSizeX, stackSizeZ, numChannels, 'uint16');
        

        for cIndex = 1:numChannels
        
            % Find the correct file for the current time point, channel, and z-slice
matchIdx = find([fileMeta.T] == (tIndex - 1) & [fileMeta.Ch] == (cIndex - 1), 1);
if isempty(matchIdx)
    error('Missing file for T=%d, Ch=%d', tIndex - 1, cIndex - 1);
end
fullFileName = fullfile(inputFolder, fileMeta(matchIdx).name);
            intermediateStack(:, :, :, cIndex) = readtiff_parallel(fullFileName);
        end

        for zIndex = 1:stackSizeZ
            for cIndex = 1:numChannels
            
                % % Find the correct file for the current time point, channel, and z-slice
                % filePattern = sprintf('*_T%04d_Ch%d.tif', tIndex - 1, cIndex-1); % Updated to match zero-padded timepoints
                % matchingFiles = dir(fullfile(inputFolder, filePattern));
                % 
                % if isempty(matchingFiles)
                %     error('Missing file for T=%04d, Ch=%d', tIndex - 1, cIndex);
                % end
                % 
                % % Read and write each slice
                % fullFileName = fullfile(inputFolder, matchingFiles(1).name);
                % slice = readtiff_parallel(fullFileName);

                t.setTag(tagstruct);
                t.write(intermediateStack(:, :, zIndex, cIndex));
                t.writeDirectory();
            end
        end
    end

    % Close the TIFF file
    t.close();
    elapsedTime = toc;
    fprintf('5D stack saved successfully with pixel spacing and frame interval in %.2f seconds.\n', elapsedTime);
end
