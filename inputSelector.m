function [inputFolder, filePaths]=inputSelector()
    inputFolder = uigetdir('', 'Select images folder');
    num_channels = inputdlg("Number of channels", "channels", 1, {'1'});
    filePaths = {};
    defaultPath=inputFolder;
    for i = 1:str2double(num_channels{1})
        text = sprintf('Select PSF file for channel %d\n',i-1);
        [file, path] = uigetfile('*.tif', text, defaultPath);
        if file == 0
            disp('User canceled file selection.');
            break;  % Exit the loop if the user cancels
        end
    defaultPath = fullfile(path, file);
    filePaths{i} = defaultPath;  % Store the full file path
    end
end
