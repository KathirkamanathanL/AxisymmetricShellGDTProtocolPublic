function convertPTS2BIN(tower)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 21.58 on 28/04/2025


% Function generates .bin files if they do not exist

% Finds all files which are .bin
files = dir(fullfile(strcat('Input_', tower),'*.bin'));

% If there are no .bin files
if isempty(files)
    
    disp('.bin files not found.')

    % Finds all files which are .pts
    files = dir(fullfile(strcat('Input_', tower),'*.pts'));

    % Looping through all .pts files
    for j =  1:length(files)
        fprintf('Binning %g of %g .pts files\n', j, length(files))
        splitted = split(files(j).name, '.');
        eval(sprintf('! %sPTS2BIN.exe %s', strcat('Common', filesep), strcat('Input_', tower, filesep, splitted{1}))) % Running PTS2BIN.exe to generate .bin file
    end

else
    disp('.bin files found. (Make sure all .bin files are generated!)')
end
