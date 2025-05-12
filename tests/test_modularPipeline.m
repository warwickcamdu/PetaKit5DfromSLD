classdef test_modularPipeline < matlab.unittest.TestCase
    methods (Test)
        function testPipelineWithSld(testCase)
            testCase.verifyPipeline('sld', {
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture2_0point6_zSpacing_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture2_0point6_zSpacing_deskew_MAX.tif'
            },{ 
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing'
                'example-data-2C-2T-3Z-1Series_Capture2_0point6_zSpacing'
            });
        end

        function testPipelineWithSldy(testCase)
            testCase.verifyPipeline('sldy', {
                'example-data-2C-2T-3Z-1Series_Capture_1_decondeskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_decondeskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_deskew_MAX.tif'
            }, { 
            'example-data-2C-2T-3Z-1Series_Capture_1'
            });
        end
    end

    methods (Access = private)
        function verifyPipeline(testCase, folderName, expected_files, subdirsToDelete)
            if nargin < 4
                subdirsToDelete = {};
            end

            thisFile = mfilename('fullpath');
            [thisDir, ~, ~] = fileparts(thisFile);
            psfPath = fullfile(thisDir, 'data', 'PSFs');
            inputPath = fullfile(thisDir, 'data', 'Inputs', folderName);
            %expectedDir = fullfile(thisDir, 'data', 'Outputs', folderName);

            % Delete output subdirectories if requested
            for i = 1:numel(subdirsToDelete)
                subdir = fullfile(inputPath, subdirsToDelete{i});
                if isfolder(subdir)
                    rmdir(subdir, 's');
                end
            end

            % Remove any old output files
            for i = 1:numel(expected_files)
                filePath = fullfile(inputPath, expected_files{i});
                if isfile(filePath)
                    delete(filePath);
                end
            end

            % Run pipeline
            modularPipeline(psfPath, inputPath);

            % Validate results
            for i = 1:numel(expected_files)
                file = expected_files{i};
                actualFile = fullfile(inputPath, file);
        
                %expectedFile = fullfile(expectedDir, file);

                testCase.assertTrue(isfile(actualFile), ...
                    sprintf('Missing output file: %s', actualFile));
                %testCase.assertTrue(isfile(expectedFile), ...
                 %   sprintf('Missing expected file: %s', expectedFile));

                % Uncomment this if pixel-level comparison is needed
                % actualImage = readtiff_parallel(actualFile);
                % expectedImage = readtiff_parallel(expectedFile);
                % testCase.assertTrue(isequal(actualImage, expectedImage), ...
                %     sprintf('Mismatch in file contents: %s', file));
            end
        end
    end
end
