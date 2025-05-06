function [THETg, Zg, RBIG] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 17.51 on 02/05/2025


% Function reconstructs shell outer surface

% Additional parameters
nRandInds = 1000; % Number of random points to consider for point spacing
searchRadius = 50; % Search radius initially used for 1NN for characteristic spacing
kNear = 2; % Minimum number of points algorithm must find to get 1NN for characteristic spacing
spacingPrctile = 95; % Percentile used to determine characteristic spacing

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute cartesian coordinates x, y of all points.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Applying coordinate transformation
XX = RR_REM.*cos(TT_REM);
YY = RR_REM.*sin(TT_REM);
ZZ = ZZ_REM;
RR = RR_REM;
TT = TT_REM;
ZZ = ZZ - botExtent;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Select nRandom points of the can to compute nearest neighbours for.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

randInds = randsample(length(TT_REM),nRandInds); % Random idexes

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for point in nRandom points do
%     Compute nearest neighbour distance from point to full can dataset and append to minDist.
% end for
% Set gridSpacing to 95th percentile of minDist.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Setting up octree
npointsbin = 1e4; % Number of points per leaf in octree
points = [XX YY ZZ;];
initialBounds = [min(min(XX)) min(min(YY)) min(min(ZZ)) max(max(XX)) max(max(YY)) max(max(ZZ)); ];
initialBounds = [min(initialBounds(:,1))-10 min(initialBounds(:,2))-10 min(initialBounds(:,3))-10 max(initialBounds(:,4))+10 max(initialBounds(:,5))+10 max(initialBounds(:,6))+10];
[binParents, binCorners, pointBins] = createOctTree(points, npointsbin, initialBounds);


searchRadiusUsed = searchRadius; % Actual searchRadius used in loop
minDist = zeros(1,nRandInds); % Array to store minimum distances
tic
for j = 1:nRandInds
    while true
        % Finding nearest neighbours
        [dist, Idx] = searchOctTreeRadius([XX(randInds(j)) YY(randInds(j)) ZZ(randInds(j))], searchRadiusUsed, points, binParents, binCorners, pointBins);

        % If at least kNear points are found
        if length(Idx) >=kNear
            [dist, sortInds] = sort(dist); Idx = Idx(sortInds); % Sorting into ascending order
            minDist(j) = dist(2); % Extracting 1NN distance
            break
        else
            % This else clause is very unlikely to be reached provided good
            % coverage with no downsampling
            searchRadiusUsed = searchRadiusUsed * 2;
            fprintf('Doubling search radius, j = %g\n', j)
        end
    end
end
toc

usedSpacing = prctile(minDist, spacingPrctile); % Setting the spacing to be used to a percentile of 1NN

%% Best Fitting Cone

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Solve for the top, rt, and bottom, rb, radius and horizontal, xbc, and vertical, ybc coordinate, 
% of centroid of best fitting truncated cone.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

xm = [0 0 0 0]; % Initial guess

options = optimoptions('lsqnonlin');
options.MaxFunctionEvaluations = 500000;
options.FunctionTolerance = 1e-14;
options.StepTolerance = 1e-14;
options.MaxIterations = 5000;
tic

f = @(b) basicCone(b, [XX YY ZZ]);
xs = lsqnonlin(f,xm,[],[], [], [],[],[],[],options);
toc

%% Defining mesh resolution

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Define vertical coordinate zng for the nominal grid.
% Estimate the variation of the radius rh with height with a linear function using rt and rb.
% Define circumferential coordinate θng for the nominal grid.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

z = linspace(0, topExtent-botExtent, ceil((topExtent-botExtent)/usedSpacing)); % Vertical point spacing

% Fitting a straight line to model the straight meridian
z1 = min(ZZ);
z2 = max(ZZ);

r1 = xs(3); % Bottom of can radius
r2 = xs(4); % Top of can radius

m = (r2-r1)/(z2-z1);
c = r1-m*z1;
rloc = m*z+c;

% Circumferential spacing - smallest radius used as grid spacing needs to
% be larger than usedSpacing
theta = linspace(0, 2*pi,ceil((2*pi)/(usedSpacing/min(rloc)))); theta = theta(1:end-1);


%% Creating nominal grid

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Create θg, Zg which are the meshgrid representations of θng and zng respectively.
% Create initial nominal Xg and Yg from θg and rh and shift to best fit centre of point cloud.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

[THETg, Zg] = meshgrid(theta, z);
RBIG = rloc'.*ones(size(THETg));
XNEW = RBIG.*cos(THETg); YNEW = RBIG.*sin(THETg); ZNEW = Zg;

% Shifting nominal grid to centre of data.
XNEW = XNEW + xs(1);
YNEW = YNEW + xs(2);

%% Gridding

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute searchRadius using gridSpacing and searchRadiusGridSpacingSF.
% 
% for circumferential coordinate in θng do
%     for vertical coordinate in zng do
%         Reset searchRadius to initial value.
%         while true do
%             Find nearest neighbours to chosen node point within the defined searchRadius.
%             if no nearest neighbours are found then
%                 Double search radius.
%             else
%                 Apply Equation (2) and compute updated grid radial value and store in Rg.
%                 break
%             end if
%         end while
%     end for
% end for
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

searchRadius = usedSpacing*searchRadiusGridSpacingSF;
kNear = 1;

% Setting up octree
points = [XX YY ZZ;];
initialBounds = [min(min(XNEW)) min(min(YNEW)) min(min(ZNEW)) max(max(XNEW)) max(max(YNEW)) max(max(ZNEW)); min(min(XX)) min(min(YY)) min(min(ZZ)) max(max(XX)) max(max(YY)) max(max(ZZ)); ];
initialBounds = [min(initialBounds(:,1))-10 min(initialBounds(:,2))-10 min(initialBounds(:,3))-10 max(initialBounds(:,4))+10 max(initialBounds(:,5))+10 max(initialBounds(:,6))+10];
[binParents, binCorners, pointBins] = createOctTree(points, npointsbin, initialBounds);

for j = 1:size(XNEW,1)
    fprintf('Gridding row %g out of %g\n', j, size(XNEW,1))
    for k = 1:size(XNEW,2)
        searchRadiusUsed = searchRadius;
        % Trialling a search radius until points are found
        while true
            [dist, Idx] = searchOctTreeRadius([XNEW(j,k) YNEW(j,k) ZNEW(j,k)], searchRadiusUsed, points, binParents, binCorners, pointBins);
            if length(Idx) >=kNear
                [dist, sortInds] = sort(dist); Idx = Idx(sortInds); % Sorting into ascending order
                RBIG(j,k) = sum((RR(Idx)./(dist.^p)))/sum(1./(dist.^p));
                break
            else
                % This is very unlikely to be reached provided good
                % coverage with no downsampling
                searchRadiusUsed = searchRadiusUsed * 2;
                fprintf('Doubling search radius, j = %g, k = %g\n', j, k)
            end
        end
    end
end

%% Plotting surface - not essential for algorithm and not in pseudocode
TT(TT<0) = TT(TT<0) + 2*pi;
figure;
hold on
scatter3(TT, ZZ, RR, 0.1, '.k')
surf(THETg, Zg, RBIG)
grid on
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
title(sprintf('Surface reconstruction of strake at %g - %g m, before smoothing', Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])

figure;
hold on
surf(THETg, Zg, RBIG)
grid on
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
title(sprintf('Surface reconstruction of strake at %g - %g m, before smoothing', Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])

%% Smoothing final surface

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute gaussian filter size using nominal can properties.
% Apply gaussian filter with standard deviation smoothingFilterStd to Rg for smoothing.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

bendingWaveLength = 2.444*sqrt(mean([R0_BOTS(S) R0_TOPS(S)])*THICKS(S)); % Bending halfwavelength
bendingBasedFilterSize = floor(bendingWaveLength/5/usedSpacing); % Assuming the maximum filter size is approximately 1/10 a bending halfwave length
bendingBasedFilterSize = 2*floor(bendingBasedFilterSize/2)+1; % Ensuring filter size is an odd number
RBIG = imgaussfilt(RBIG,smoothingFilterStd,'FilterSize',max([bendingBasedFilterSize 3])); % Appling filter with minimum size of 3.

%% Plotting filtered surface after smoothing - not essential for algorithm and not in pseudocode
TT(TT<0) = TT(TT<0) + 2*pi;
figure;
hold on
scatter3(TT, ZZ, RR, 0.1, '.k')
surf(THETg, Zg, RBIG)
grid on
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
title(sprintf('Surface reconstruction of strake at %g - %g m, after smoothing', Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])

figure;
hold on
surf(THETg, Zg, RBIG)
grid on
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
title(sprintf('Surface reconstruction of strake at %g - %g m, after smoothing', Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])


end

function [perpImp] = basicCone(b, x)

% b(1) = X1
% b(2) = Y1
% b(3) = Botradius
% b(4) = Topradius

X = x(:,1); % Raw X coords
Y = x(:,2); % Raw Y coords
Z = x(:,3); % Raw Z coords
Z = Z - min(Z);

H = max(Z);

coords = [X-b(1) Y-b(2) Z]; % Shifting origin


R = (b(3) - b(4))/(H - 0)*(H - Z)+b(4);
beta = atan2((b(3) - b(4)), H);

radialImp = R - sqrt((coords(:,1)).^2 + (coords(:,2)).^2);
perpImp = abs(radialImp.*cos(beta));

end