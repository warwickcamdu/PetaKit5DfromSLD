% File: testmodularPipeline.m
classdef test_modularPipeline < matlab.unittest.TestCase
    methods (Test)
        function testPipelineWithSld(testCase)
            thisFile = mfilename('fullpath');
            [thisDir, ~, ~] = fileparts(thisFile);
            psfPath = fullfile(thisDir,'data', 'PSFs');
            inputPath = fullfile(thisDir, 'data', 'Inputs', 'sld');
            % Call the function being tested
            modularPipeline(psfPath, inputPath); 
            expected_files = {
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_decondeskew_MAX.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew.tif'
                'example-data-2C-2T-3Z-1Series_Capture1_0point5_zSpacing_deskew_MAX.tif'
               };
            % Assert that the result is correct
            %testCase.verifyEqual(result, 5);
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
