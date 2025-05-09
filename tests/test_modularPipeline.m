% File: testmodularPipeline.m
classdef test_modularPipeline < matlab.unittest.TestCase
    methods (Test)
        function testPipelineWithSld(testCase)
            thisFile = mfilename('fullpath');
            [thisDir, ~, ~] = fileparts(thisFile);
            psfPath = fullfile(thisDir,'data', 'PSFs');
            inputPath = fullfile(thisDir, 'data', 'Inputs', 'sld');
            expected_files = {
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point6_zSpacing_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point6_zSpacing_deskew_MAX.tif'
               };
            expectedDir=fullfile(thisDir, 'data', 'Outputs', 'sld');
            tifDir = fullfile(inputPath,'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing');
            if isfolder(tifDir)
                                   rmdir(tifDir,'s');
            end
              tifDir = fullfile(inputPath,'example-data-2C-2T-3Z-1Series_Capture1_0point6_zSpacing');
            if isfolder(tifDir)
                rmdir(tifDir,'s');
            end
                % Pre-test cleanup: delete expected files if they already exist
    for i = 1:numel(expected_files)
        filePath = fullfile(inputPath, expected_files{i});
        if isfile(filePath)
            delete(filePath);
        end
    end
       
            % Call the function being tested
            modularPipeline(psfPath, inputPath); 
            
            % Assert that the result is correct
                for i = 1:numel(expected_files)
        file = expected_files{i};
        actualFile = fullfile(inputPath, file);
        expectedFile = fullfile(expectedDir, file);

        testCase.assertTrue(isfile(actualFile), ...
            sprintf('Missing output file: %s', actualFile));
        testCase.assertTrue(isfile(expectedFile), ...
            sprintf('Missing expected file: %s', expectedFile));

        actualImage = readtiff_parallel(actualFile);
        expectedImage = readtiff_parallel(expectedFile);

        testCase.assertTrue(isequal(actualImage, expectedImage), ...
            sprintf('Mismatch in file contents: %s', file));
    end
        end
        
        function testPipelineWithSldy(testCase)
            thisFile = mfilename('fullpath');
            [thisDir, ~, ~] = fileparts(thisFile);
            psfPath = fullfile(thisDir,'data', 'PSFs');
            inputPath = fullfile(thisDir, 'data', 'Inputs', 'sldy');
            % Simulate calling the function using those folders
            modularPipeline(psfPath, inputPath); % assuming you modify the function to accept inputs
            expected_files = {
                'example-data-2C-2T-3Z-1Series_Capture_1_decondeskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_decondeskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture_1_deskew_MAX.tif'
               };
            % Add assertions here
        end
    end
end
