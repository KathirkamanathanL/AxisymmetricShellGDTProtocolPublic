function [regParams, bestFitParamsPreReg, iUse, nClouds] = S1_EnhancedRegistration(minZ, maxZ, units, nperc, tower, numIter, cutoffReg)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 14.52 on 06/05/2025


% Function performs enhanced registration

%% Loading and plotting data

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute number of data clouds, numClouds.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

files = dir(fullfile(strcat('Input_', tower),'*.bin')); % Finds all files which are .bin
nClouds = length(files); % Number of clouds
COBJS = []; % Array to store S1_CLOUD objects

% Finding minimum shell radius
run(fullfile(strcat('Input_', tower),strcat('S0_Data_', tower)))
rmin = min([R0_BOTS R0_TOPS])/1e3; % Minimum radius of shell

% Loading data, initial preprocessing, object creation
for j = 1:nClouds % Looping through all files
    if j == 1
        COBJS = S1_CLOUD(files(j).name, j); % Creating a S1_CLOUD object and cloud array
    else
        COBJS(end+1) = S1_CLOUD(files(j).name, j); % Adding to the cloud array
    end

    COBJS(j) = COBJS(j).CLOUDDOWNRND(nperc); % Downsampling clouds
    COBJS(j) = COBJS(j).CLOUDZLIM(minZ,maxZ); % Limiting the range of the cloud
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
end

% Plotting the data
figure;
hold on
for j = 1:nClouds % Looping through all clouds
    scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
end
xlabel('$\theta$ [rad]', 'interpreter', 'latex')
ylabel(['$z$ [', units ,']'], 'interpreter', 'latex')
zlabel(['$\rho$ [', units ,']'], 'interpreter', 'latex')
title('Cylindrical coordinates WTST', 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
grid on

figure;
hold on
for j = 1:nClouds % Looping through all clouds
    scatter3(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, 1, 'marker', '.')
end
xlabel(['$x$ [', units ,']'], 'interpreter', 'latex')
ylabel(['$y$ [', units ,']'], 'interpreter', 'latex')
zlabel(['$z$ [', units ,']'], 'interpreter', 'latex')
title('Cartesian coordinates WTST', 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
grid on

%% Applying best fit cone correction

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Fit a best-fitting truncated cone to all clouds at once, apply transformation and store 
% parameters in bestFitParamsPreReg.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Isolating desired region
X = []; Y = []; Z = [];

for j = 1:nClouds % Looping through all files
    X(end+1:end+length(COBJS(j).X)) = COBJS(j).X;
    Y(end+1:end+length(COBJS(j).Y)) = COBJS(j).Y;
    Z(end+1:end+length(COBJS(j).Z)) = COBJS(j).Z;
end

sel = Z < minZ | Z > maxZ;
X(sel) = []; Y(sel) = []; Z(sel) = [];

% Computing best fit cone parameters for entire shell
[bestFitParamsPreReg] = bestFitCone(X, Y, Z, minZ, tower, 1e6);

% Applying transformation to each scan
for j = 1:nClouds % Looping through all clouds
    [COBJS(j).X, COBJS(j).Y, COBJS(j).Z, COBJS(j).T, COBJS(j).R] = bestFitConeTransform(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, minZ, bestFitParamsPreReg);
end

%% Iterative loop for ICP

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise regParams to store registration parameters.
% for iteration in (1 to numIter + 1) do
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Initialising arrays
registrationSolutions = []; % Array to store registration constants from each iteration
ptp99 = zeros(numIter+1,1); % Array to store point to plane metric

for j = 1:nClouds
    COBJS(j).XNEW = COBJS(j).X;
    COBJS(j).YNEW = COBJS(j).Y;
    COBJS(j).ZNEW = COBJS(j).Z;
end

for iterCounter=1:numIter+1
    fprintf('Iteration = %g\n', iterCounter)
    for j = 1:nClouds
        COBJS(j).X = reshape(COBJS(j).XNEW, length(COBJS(j).XNEW),1);
        COBJS(j).Y = reshape(COBJS(j).YNEW, length(COBJS(j).XNEW),1);
        COBJS(j).Z = reshape(COBJS(j).ZNEW, length(COBJS(j).XNEW),1);
    end


    %% Finding extent of each cloud

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % 	for cloud in clouds do
    % 		Divide cloud into circumferential bins and identify if there are points in each bin.
    %       Extend bins circumferentially by necessary amount.
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    CLOUDExtent = {}; % Stores the bins the cloud is present in
    cloudConn = []; % Connections between clouds considered i.e. [1 2; 1 3....]
    CONNECTIONS = {}; % All connections between clouds
    CONNECTIONSdist = {};  % Distance between all connections


    thetOvBins = linspace(0, 2*pi, 101); % Points divinding circumference to determine extent of each scan

    for j = 1:nClouds

        thetExtend = cutoffReg/rmin; % Maximum angle around shell that cutoff distance corresponds to
        thetOvBinSize = thetOvBins(2) - thetOvBins(1); % Size of each theta bin when looking for circumferential extent

        numExtend = ceil(thetExtend/thetOvBinSize); % Number bins to extend circumferential extent accounting for cutoff size
        extentCloud = zeros(length(thetOvBins)-1, 1); % Boolean array storing whether a bin has points in it or not

        % Cylindrical coordinates
        T = atan2(COBJS(j).Y, COBJS(j).X);
        R = sqrt(COBJS(j).X.^2 + COBJS(j).Y.^2);
        T(T < 0) = T(T < 0) + 2*pi;

        % Looping through bins and storing whether there are any points in extentCloud
        for k = 1:length(thetOvBins)-1

            if k == 1
                sel = T >= thetOvBins(k) & T <= thetOvBins(k+1);
            else
                sel = T > thetOvBins(k) & T <= thetOvBins(k+1);
            end

            if any(sel)
                extentCloud(k) = 1;
            end

        end

        % Extending extentCloud variable to account for cutoffReg for acceptable pairings in ICP
        extentInds = find(extentCloud); % Index of bins with points in them
        extendedInds = extentInds; % Index of bins after extension by numExtend
        for j = 1:numExtend
            extendedInds = union(extendedInds,extentInds-j);
            extendedInds = union(extendedInds,extentInds+j);
        end

        % Ensuring the bins after extending are within appropriate ranges
        extendedInds(extendedInds <= 0) = extendedInds(extendedInds <= 0) + length(extentCloud);
        extendedInds(extendedInds > length(extentCloud)) = extendedInds(extendedInds > length(extentCloud)) - length(extentCloud);
        extendedInds = unique(extendedInds);

        % Storing the extended extent of each cloud
        CLOUDExtent{end+1} = extendedInds;
    end


    %% Finding all connections to be used in registration

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   for firstCloudIndex in (1 to numClouds) do
    % 		for secondCloudIndex in (firstCloudIndex + 1 to numClouds) do
    % 			for overlappingBins between firstCloudIndex and secondCloudIndex clouds do
    % 				for point in overlappingBin from cloud with secondCloudIndex do
    % 					Find nearest neighbour of point, to points in cloud with firstCloudIndex
    % 				end for
    % 			end for
    %           Remove nearest neighbours with distances more than cutoffReg.
    %           Store correspondences in cell array called connections.
    % 		end for
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    knear = 1; % Number of desired nearest neighbours
    npointsbin = ceil(length(COBJS(1).X)/1e2); % Number of points in octree bin

    tic
    for j = 1:nClouds
        % Creating octree
        points = [COBJS(j).X COBJS(j).Y COBJS(j).Z];
        initialBounds = [min(X)-1 min(Y)-1 min(Z)-1 max(X)+1 max(Y)+1 max(Z)+1];
        [binParents, binCorners, pointBins] = createOctTree([COBJS(j).X COBJS(j).Y COBJS(j).Z], npointsbin, initialBounds);

        for m = j+1:nClouds
            fprintf('1NN for clouds %g and %g \n', j, m)
            cloudConn(end+1,1:2) = [j m]; % Storing cloud pair being considered

            % Finding overlaps between two clouds
            intersectInds = intersect(CLOUDExtent{j}, CLOUDExtent{m});

            % Radial coordinate of cloud m
            T = atan2(COBJS(m).Y, COBJS(m).X);
            T(T < 0) = T(T < 0) + 2*pi;

            % Initialising arrays
            dists = []; nearestInds = []; allPoints = [];

            % Looping through overlapping bins
            for n = 1:length(intersectInds)

                % Identifying points in bin
                if intersectInds(n) == 1
                    sel = T >= thetOvBins(intersectInds(n)) & T <= thetOvBins(intersectInds(n)+1);
                else
                    sel = T > thetOvBins(intersectInds(n)) & T <= thetOvBins(intersectInds(n)+1);
                end
                searchOverlapInds = find(sel);
                allPoints(end+1:end+length(searchOverlapInds)) = searchOverlapInds;

                % Initialising arrays for parfor loop
                distsSubArray = zeros(length(searchOverlapInds),1);
                nearestIndsSubArray = distsSubArray;

                % Performing nearest neighbour search
                parfor k = 1:length(searchOverlapInds)
                    [distsSubArray(k), nearestIndsSubArray(k)] = searchOctTreeKNN([COBJS(m).X(searchOverlapInds(k)) COBJS(m).Y(searchOverlapInds(k)) COBJS(m).Z(searchOverlapInds(k))], knear, points, binParents, binCorners, pointBins);
                end

                % Adding results of nearest neighbour search to arrays
                dists = [dists; distsSubArray]; nearestInds = [nearestInds; nearestIndsSubArray];
            end

            % Reshaping arrays
            dists = reshape(dists,length(dists), 1);
            nearestInds = reshape(nearestInds,length(dists),1);
            
            % Removing nearest neighbour connections outside threshold
            sel = dists > cutoffReg;
            nearestInds(sel) = [];
            allPoints(sel) = [];
            dists(sel) = [];

            % Storing feasible connections and distances between points
            CONNECTIONS{end+1} = [nearestInds allPoints'];
            CONNECTIONSdist{end+1} = dists;
        end
    end
    toc

    fprintf('Number of connections = %g\n', sum(cellfun(@numel, CONNECTIONSdist)))


    %% Estimating point to plane before registration

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   Compute point-to-plane metric for all correspondences and store 99th percentile.
    %   if iteration equals (numIter+1) then
    % 		break
    %   end if
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % Computing normals of clouds
    for j = 1:nClouds
        ptCloud = pointCloud([COBJS(j).X COBJS(j).Y COBJS(j).Z]);
        COBJS(j).normals = pcnormals(ptCloud);
    end

    % Array to store point to plane metric
    mindists = [];

    % Computing point to plane metric for two clouds
    for j = 1:size(cloudConn,1)
        c1 = cloudConn(j,:);
        c2 = c1(2); c1 = c1(1); % Relavant clouds
        conn = CONNECTIONS{j}; % Connections between clouds
        if any(conn)
            c1conn = conn(:,1); c2conn = conn(:,2);
            point2Planedists = (dot(COBJS(c1).normals(c1conn,:), [COBJS(c1).X(c1conn)-COBJS(c2).X(c2conn) COBJS(c1).Y(c1conn)-COBJS(c2).Y(c2conn) COBJS(c1).Z(c1conn) - COBJS(c2).Z(c2conn)], 2)).^2; % Point to plane
            mindists(end+1:end+length(point2Planedists)) = point2Planedists;
        end
    end

    fprintf('Point to plane metric p99 = %g\n', prctile(mindists,99))
    ptp99(iterCounter) = prctile(mindists,99);

    % If statement is neccessary as the point to plane metric is computed
    % for iteration numIter on loop counter numIter+1
    if iterCounter == numIter+1
        break
    end

    %% Solving for registration constants

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   Solve ICP minimisation problem with point-to-plane metric using lsqnonlin for all clouds 
    %   simultaneously and obtain transformation parameters.
    %   Apply solved registration transformation to clouds and append to regParams.
    % end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    nopt = 6;

    % lsqnonlin solver settings
    options = optimoptions('lsqnonlin');
    options.FunctionTolerance = 1e-6;
    options.OptimalityTolerance = 1e-6;
    options.StepTolerance = 1e-6;
    options.MaxFunctionEvaluations = 1e4;
    options.MaxIterations = 1e4;
    options.UseParallel = true;

    xm = zeros((nClouds-1)*nopt,1);
    f = @(b)ICPFunction(b, COBJS, nClouds, cloudConn, CONNECTIONS);

    % Appling lsqnonlin solver
    tic
    x = lsqnonlin(f,xm,[],[],options);
    toc

    % Adding transformations to first one
    regResult = [0 0 0 0 0 0 x']';

    % Applying registration
    for j = 1:nClouds % Looping through all clouds
        COBJS(j) = COBJS(j).CLOUDMERGE(regResult); % Applying registration
        COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).XNEW, COBJS(j).YNEW); % Computing cylindrical coordinates
        COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
    end

    % Storing registration result
    registrationSolutions = [registrationSolutions x];

end

%% Plotting point to plane metric

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Identify iUse, the iteration with smallest 99th percentile of point-to-plane metric.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

figure
plot(0:length(ptp99)-1, ptp99, '-x')
xlabel('Iteration [-]', 'Interpreter','latex')
ylabel('Point to plane metric [m]', 'Interpreter','latex')
title('Point to plane metric variation with iteration', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
grid on

% Finding iteration with the smallest point to plane metric
[~, iUse] = min(ptp99);
iUse = iUse-1;
fprintf('Most optimum fit is iteration %g\n', iUse)

% Adding zero motion for first cloud;
regParams = [zeros(6,size(registrationSolutions,2)); registrationSolutions];

%% Plotting results
% Loading data again
nperc = 0.02;
files = dir(fullfile(strcat('Input_', tower),'*.bin')); % Finds all files which are .bin
COBJS = []; % Array to store cloud objects

% Loading data, initial preprocessing, object creation
for j = 1:nClouds % Looping through all clouds
    if j == 1
        COBJS = S1_CLOUD(files(j).name, j); % Creating a S1_CLOUD object and cloud array
    else
        COBJS(end+1) = S1_CLOUD(files(j).name, j); % Adding to the cloud array
    end

    COBJS(j) = COBJS(j).CLOUDDOWNRND(nperc); % Downsampling clouds
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
end

% Applying original best fit transformation to each scan
for j = 1:nClouds % Looping through all clouds
    [COBJS(j).X, COBJS(j).Y, COBJS(j).Z, COBJS(j).T, COBJS(j).R] = bestFitConeTransform(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, minZ, bestFitParamsPreReg);
end

% Plotting initial geometry
figure
hold on
for j = 1:nClouds % Looping through all clouds
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
    COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
    scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
end
view(0,0)
grid on
zlim([3.72 3.77])
ylim([29.0 29.1])
xlabel('$\theta$ [rad]', 'interpreter', 'latex')
ylabel('$z$ [m]', 'interpreter', 'latex')
zlabel('$\rho$ [m]', 'interpreter', 'latex')
title(sprintf('Initial geometry @ %g-%g m', gca().YLim(1), gca().YLim(2)), 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')

% Plotting registered geometries
for solNum = 1:numIter
    figure
    hold on
    % Applying registration
    for j = 1:nClouds % Looping through all clouds
        COBJS(j) = COBJS(j).CLOUDMERGEITER(regParams, solNum); % Applying registration
        COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).XNEW, COBJS(j).YNEW); % Computes cylindrical coordinates
        COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
        scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
    end
    view(0,0)
    grid on
    zlim([3.72 3.77])
    ylim([29.0 29.1])
    xlabel('$\theta$ [rad]', 'interpreter', 'latex')
    ylabel('$z$ [m]', 'interpreter', 'latex')
    zlabel('$\rho$ [m]', 'interpreter', 'latex')
    title(sprintf('Iteration %g @ %g-%g m', solNum, gca().YLim(1), gca().YLim(2)), 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
end

end

%% Functions
function [mindists] = ICPFunction(b, COBJS, nClouds, cloudConn, CONNECTIONS)

% Function returns quantity to be minimsed in ICP

nopt = 6;

% Applying transformations
for j = 1:nClouds
    if j==1
        Ci = zeros(1,nopt);
    else
        Ci = b(nopt*COBJS(j-1).ID-(nopt-1):nopt*COBJS(j-1).ID);
    end

    Rxi = [1 0 0; 0 cos(Ci(1)) -sin(Ci(1)); 0 sin(Ci(1)) cos(Ci(1))]; % Rotation about x axis
    Ryi = [cos(Ci(2)) 0 sin(Ci(2)); 0 1 0; -sin(Ci(2)) 0 cos(Ci(2))]; % Rotation about Y axis
    Rzi = [cos(Ci(3)) -sin(Ci(3)) 0; sin(Ci(3)) cos(Ci(3)) 0; 0 0 1]; % Rotation about Y axis

    newcoords = Rzi*Ryi*Rxi*[COBJS(j).X COBJS(j).Y COBJS(j).Z]' + [Ci(4) Ci(5) Ci(6)]';
    newcoords = newcoords';
    COBJS(j).X = newcoords(:,1);
    COBJS(j).Y = newcoords(:,2);
    COBJS(j).Z = newcoords(:,3);

end

% Computing point to plane metric
mindists = [];
for j = 1:size(cloudConn,1)
    c1 = cloudConn(j,:);
    c2 = c1(2); c1 = c1(1); % Relavant clouds
    conn = CONNECTIONS{j}; % Connections between clouds
    if any(conn)
        c1conn = conn(:,1); c2conn = conn(:,2);
        dists = (dot(COBJS(c1).normals(c1conn,:), [COBJS(c1).X(c1conn)-COBJS(c2).X(c2conn) COBJS(c1).Y(c1conn)-COBJS(c2).Y(c2conn) COBJS(c1).Z(c1conn) - COBJS(c2).Z(c2conn)], 2)).^2;
        mindists(end+1:end+length(dists)) = dists;
    end
end
end

