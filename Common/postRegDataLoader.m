function [X, Y, Z, T, R] = postRegDataLoader(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 15.42 on 29/04/2025


% Function loads data after Stage 1 registration has been performed

files = dir(fullfile(strcat('Input_', tower),'*.bin')); % Finds all files which are .bin
nclouds = length(files); % Number of clouds
COBJS = []; % Array to store cloud objects

% Loading data, initial preprocessing, object creation
for j = 1:nclouds % Looping through all files
    if j == 1
        COBJS = S1_CLOUD(files(j).name, j); % Creating a CLOUD object and cloud array
    else
        COBJS(end+1) = S1_CLOUD(files(j).name, j); % Adding to the cloud array
    end

    COBJS(j) = COBJS(j).CLOUDDOWNRND(nperc); % Downsampling clouds
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes polar coordinate
end


% Applying original best fit transformation to each scan
for j = 1:nclouds % Looping through all clouds
    [COBJS(j).X, COBJS(j).Y, COBJS(j).Z, COBJS(j).T, COBJS(j).R] = bestFitConeTransform(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, minZ, bestFitParamsPreReg);
end

% Applying registration
for j = 1:nclouds % Looping through all clouds
    COBJS(j) = COBJS(j).CLOUDMERGEITER(regParams, iUse); % Applying registration
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).XNEW, COBJS(j).YNEW); % Computes cylindrical coordinate
end

% Isolating desired region
X = []; Y = []; Z = [];

for j = 1:nclouds % Looping through all clouds
    X(end+1:end+length(COBJS(j).XNEW)) = COBJS(j).XNEW;
    Y(end+1:end+length(COBJS(j).XNEW)) = COBJS(j).YNEW;
    Z(end+1:end+length(COBJS(j).XNEW)) = COBJS(j).ZNEW;
end

sel = Z < minZ | Z > maxZ;
X(sel) = [];Y(sel) = []; Z(sel) = [];

% Converting to cylindrical coordinates
T = atan2(Y, X);
R = sqrt(Y.^2 + X.^2);
T(T < 0) = T(T < 0) + 2*pi;

end