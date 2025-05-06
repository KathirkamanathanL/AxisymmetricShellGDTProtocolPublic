% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 13.45 on 03/05/2025


% Run the code to go through each stage of the protocol. Once a stage has
% been completed, the relevant data is saved to a .mat file in the same
% directory as the point clouds. The protocol therfore does not need to be
% run from the start each time for a dataset.

clear;
clc;
close all;

tower = 'ReleaseTower'; % Tower to be processed
startFromStage = 1; % Stage to start protocol from - 1 to 7. (3, 4 and 5 run within the same loop). Use 7 to skip to final protocol output.

%% Initial set up before registration
minZ = 8; maxZ = 53; % Range of tower being considered between minZ and maxZ
units  = 'm'; % Units to be used in plots

%% Adding folders to path
addpath(genpath('Common'))
addpath(genpath(strcat('Input_', tower)))
addpath(genpath('Stage_1'))
addpath(genpath('Stage_2'))
addpath(genpath('Stage_3'))
addpath(genpath('Stage_4'))
addpath(genpath('Stage_5'))
addpath(genpath('Stage_6'))

%% Create .bin files
convertPTS2BIN(tower) % Creates .bin files if no .bin files exist in dataset directory

%% STAGE 1: ENHANCED REGISTRATION
% #### Hyperparameters #### %
nperc = 0.01; % Remaining fraction after downsampling
numIter = 5; % Number of iterations
cutoffReg = 0.02; % Cutoff for nearest neighbour [m]

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
if startFromStage > 1 && isfile(fullfile(strcat('Input_', tower),strcat('S1_Registration_', tower, '.mat')))
    disp('STAGE 1: ENHANCED REGISTRATION - Loading')
    load(fullfile(strcat('Input_', tower),strcat('S1_Registration_', tower, '.mat')))
else
    disp('STAGE 1: ENHANCED REGISTRATION - Start')
    [regParams, bestFitParamsPreReg, iUse, nClouds] = S1_EnhancedRegistration(minZ, maxZ, units, nperc, tower, numIter, cutoffReg); % Doing registration
    save(fullfile(strcat('Input_', tower),strcat('S1_Registration_', tower, '.mat')), 'nClouds', 'regParams', 'bestFitParamsPreReg', 'iUse');
    disp('STAGE 1: ENHANCED REGISTRATION - Complete')
end

%% STAGE 2: CAN SEGMENTATION
% #### Hyperparameters - Global radial profile generation #### %
nperc = 0.1;
expectedBeadWidth = 0.08; % m - Setting expected bead width
numCircWindows = 100; % Number of circumferential windows
subIntervals = 10; % Number of subIntervals. This controls the window size. Must be at least 2.
% iUse = 2; % Computed with smallest 99th percentile point to plane metric. Uncomment to overule and manually choose a different iteration

% #### Hyperparameters - Joint ROI identification (Some hyperparameters carry over from global radial profile) #### %
nJoints = 16; % Number of joints in data

% #### Hyperparameters - Joint segmentation #### %
flangeTransitionJointIndex = [9]; % Indexes of jointCoordinates where flange-flange joint exists
flangePadding = 10; % [mm] Amount to extend region identified as flange-flange joint
weldPadding = 5; % [mm] Amount to extend region identified as welded joint

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
if startFromStage > 2 && isfile(fullfile(strcat('Input_', tower),strcat('S2_CanSegmentation_', tower, '.mat')))
    disp('STAGE 2: CAN SEGMENTATION - Loading')
    load(fullfile(strcat('Input_', tower),strcat('S2_CanSegmentation_', tower, '.mat')))
else

    disp('STAGE 2: CAN SEGMENTATION - Start')

    % Global radial profile generation 
    disp('STAGE 2: Global radial profile generation ')
    [ZMEDIANS, RMEDIANS] = S2_GlobalRadialProfileGeneration(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower, expectedBeadWidth, numCircWindows, subIntervals);
    
    % Joint ROI identification
    disp('STAGE 2: Joint ROI identification')
    [jointCoordinates] = S2_JointROIIdentification(ZMEDIANS, RMEDIANS, nJoints, expectedBeadWidth);

    % Joint segmentation
    disp('STAGE 2: Joint segmentation')
    [JOINTWIDTH, JOINTBOUNDS] = S2_JointSegmentation(ZMEDIANS, RMEDIANS, jointCoordinates, flangeTransitionJointIndex, flangePadding, weldPadding);
    save(fullfile(strcat('Input_', tower),strcat('S2_CanSegmentation_', tower, '.mat')), 'ZMEDIANS', 'RMEDIANS', 'jointCoordinates', 'JOINTWIDTH', 'JOINTBOUNDS');

    disp('STAGE 2: CAN SEGMENTATION - Complete')
end

%% Preparing data for looping through cans
cloudStrakes = 1:15; % Strakes to be meshed from cloud data - Check if within minZ and maxZ
nperc = 0.02; % Fraction of data to use for remaining stages

% If a later start stage is requested and the relevant .mat file exists, no need to load point cloud
if startFromStage > 5 && isfile(fullfile(strcat('Input_', tower),strcat('S5_Surface_', tower, '.mat')))
    disp('')
else
    
    % Loading data
    [X, Y, Z, ~, ~] = postRegDataLoader(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower);
    
    % Performing best fit cone
    [bestFitParams] = bestFitCone(X, Y, Z, minZ, tower, 1e6);
    [X, Y, Z, T, R] = bestFitConeTransform(X, Y, Z, minZ, bestFitParams);
    
    % Loading tower details
    run(fullfile(strcat('Input_', tower),strcat('S0_Data_', tower)))
    
    sel = Z_BOTS/1e3 > minZ & Z_BOTS/1e3 < maxZ;
    firstStrakeInData = find(sel,1); % ID of first strake in data i.e bottom of which strake is first found in jointCoordinates
    JOINTBOUNDS = [zeros(firstStrakeInData-1,2); JOINTBOUNDS]; % Adding extra numbers to JOINTBOUNDS to aid indexing and ensure the nth joint corresponds to the nth can if a different region of shell is considered via different minZ and maxZ
    
    % Applying unit conversion of point cloud to mm
    X = X*1000; Y = Y*1000; Z = Z*1000; R = R*1000;
    JOINTBOUNDS = JOINTBOUNDS*1e3;
    N1 = zeros(size(Z_BOTS));
    N2 = N1;
end

%% Looping through cans and performing stages 3, 4 and 5 for each can

% #### Hyperparameters - STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL #### %
cutOff = 1; % [mm] Cutoff distance to determine whether a point is an outlier
windowWidthSF = 3; % Scale factor applied to expectedBeadWidth to determine window width

% #### Hyperparameters - STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION #### %
p = 1; % Power used in scatter points method
searchRadiusGridSpacingSF = 4; % Scale factor applied to grid spacing to determine search radius for scatter points method
smoothingFilterStd = 10; % Standard deviation amplitude for smoothing filter

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stages
if startFromStage > 5 && isfile(fullfile(strcat('Input_', tower),strcat('S5_Surface_', tower, '.mat')))
    disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Loading')
    disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Loading')
    disp('STAGE 5: SURFACE INPAINTING - Loading')
    load(fullfile(strcat('Input_', tower),strcat('S5_Surface_', tower, '.mat')))
else
    for S = 1:length(Z_BOTS)
        disp(['Creating Strake ', num2str(S)])
    
        if ismember(S, cloudStrakes)
            TT = T;
            RR = R;
            ZZ = Z;
            XX = X;
            YY = Y;
           
            botExtent = JOINTBOUNDS(S,2); topExtent = JOINTBOUNDS(S+1,1); % Top and bottom of can
            JOINTCUTOFFTOP = ((JOINTBOUNDS(S+1,1)+JOINTBOUNDS(S+1,2))/2-JOINTBOUNDS(S+1,1)); % Distance from top of can to middle of above joint
            JOINTCUTOFFBOT = (JOINTBOUNDS(S,2)-((JOINTBOUNDS(S,1)+JOINTBOUNDS(S,2))/2)); % Distance from bottom of can to middle of below joint
    
            % Identifying can region
            sel = Z < botExtent | Z > topExtent;
            TT(sel) = [];
            ZZ(sel) = [];
            RR(sel) = [];
            XX(sel) = [];
            YY(sel) = [];
    
            % STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL
            disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Start')
            tbinSize = expectedBeadWidth*windowWidthSF*1e3; zbinSize = tbinSize; % Window sizes for outlier removal
            if topExtent-botExtent < zbinSize; zbinSize = topExtent-botExtent-5; end % Accounting for if can is smaller than zbinSize
            [TT_REM, ZZ_REM, RR_REM, ~, ~, ~, ~, ~, ~, ~] = S3_SystematicAndRandomOutlierRemoval(TT, ZZ, RR, (R0_BOTS(S) + R0_TOPS(S))/2, tbinSize, zbinSize, cutOff);
            title(sprintf('Outlier removal of strake at %g - %g m',Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
            disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Complete')

            % STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION
            disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Start')
            [THETg, Zg, RBIG] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd);
            disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Complete')

            % STAGE 5: SURFACE INPAINTING
            disp('STAGE 5: SURFACE INPAINTING - Start')
            [THETg, Zg, RBIG] = S5_SurfaceInpainting(THETg, Zg, RBIG, JOINTCUTOFFBOT, JOINTCUTOFFTOP);
            disp('STAGE 5: SURFACE INPAINTING - Complete')
    
            figure;
            hold on
            surf(THETg, Zg, RBIG)
            grid on
            xlabel('$\theta$ [rad]','interpreter','latex')
            ylabel('$z$ [mm]','interpreter','latex')
            zlabel('$\rho$ [mm]','interpreter','latex')
            title(sprintf('Surface inpainting of strake at %g - %g m',Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
            set(gca,'TickLabelInterpreter','latex')
            set(gcf,'Position',[297 275 863 516])
            set(gca, 'View',[-16.7552 56.8675])
    
            % Computing cartesian coordinates of can
            X3D = RBIG.*cos(THETg); Y3D = RBIG.*sin(THETg); Z3D = Zg;
            Z3D = Z3D + botExtent-JOINTCUTOFFBOT;
    
            n1 = size(X3D,2); % Number of nodes circumferentially
            n2 = size(X3D,1); % Number of nodes vertically
            N1(S) = n1; N2(S) = n2;
    
            % Storing reconstruction outputs
            CANS.(sprintf('S%g',S)) = [];
    
            CANS.(sprintf('S%g',S)).THETg = THETg;
            CANS.(sprintf('S%g',S)).Zg = Zg;
            CANS.(sprintf('S%g',S)).RBIG = RBIG;
    
            CANS.(sprintf('S%g',S)).X3D = X3D;
            CANS.(sprintf('S%g',S)).Y3D = Y3D;
            CANS.(sprintf('S%g',S)).Z3D = Z3D;
    
            CANS.(sprintf('S%g',S)).N1 = n1;
            CANS.(sprintf('S%g',S)).N2 = n2;
    
        else
            N1(S) = 0; N2(S) = 0;
    
            CANS.(sprintf('S%g',S)) = [];
    
            CANS.(sprintf('S%g',S)).THETg = [];
            CANS.(sprintf('S%g',S)).Zg = [];
            CANS.(sprintf('S%g',S)).RBIG = [];
    
            CANS.(sprintf('S%g',S)).X3D = [];
            CANS.(sprintf('S%g',S)).Y3D = [];
            CANS.(sprintf('S%g',S)).Z3D = [];
    
            CANS.(sprintf('S%g',S)).N1 = 0;
            CANS.(sprintf('S%g',S)).N2 = 0;
        end
    
        disp(' ')
        
    end
    
    % Plotting all reconstructed cans on one figure
    figure
    hold on
    for S = 1:length(Z_BOTS)
        surf(CANS.(sprintf('S%g',S)).THETg, CANS.(sprintf('S%g',S)).Z3D, CANS.(sprintf('S%g',S)).RBIG)
    end
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    grid on
    set(gca,'TickLabelInterpreter','latex')
    set(gcf,'Position',[297 275 863 516])
    set(gca, 'View',[-16.7552 56.8675])

    save(fullfile(strcat('Input_', tower),strcat('S5_Surface_', tower, '.mat')), 'cloudStrakes', 'CANS', 'N1', 'N2', 'Z_BOTS', 'Z_TOPS', 'R0_BOTS', 'R0_TOPS', 'THICKS')
end

%% STAGE 6: PROJECTION TO ARBITRARY MESHES
% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
if startFromStage > 6 && isfile(fullfile(strcat('Input_', tower),strcat('S6_Mesh_', tower, '.mat')))
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Loading')
    load(fullfile(strcat('Input_', tower),strcat('S6_Mesh_', tower, '.mat')))
else
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Start')
    [X_MESH, Y_MESH, Z_MESH, THET_MESH, R_MESH, N1_MESH, N2_MESH] = S6_ProjectionToArbitraryMeshes(CANS, cloudStrakes, Z_BOTS, Z_TOPS, R0_BOTS, R0_TOPS, THICKS);
    save(fullfile(strcat('Input_', tower),strcat('S6_Mesh_', tower, '.mat')), 'cloudStrakes', 'X_MESH', 'Y_MESH', 'Z_MESH', 'THET_MESH', 'R_MESH', 'N1_MESH', 'N2_MESH', 'Z_BOTS', 'Z_TOPS', 'R0_BOTS', 'R0_TOPS', 'THICKS')
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Complete')
end
