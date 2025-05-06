% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 20.43 on 06/05/2025


% Code measures the circumferential dimple tolerance. 

% Measurements are to be taken with the gauge placed on the outside of the
% shell only. There is an option to use the gauge on the inside of the
% shell, but is not accurate as the gauge will always pass through the
% shell wall in such situations.

% An additional check is provided via the selfIntersectionCheck function
% which makes sure no part of the curved gauge (especially from gaugeStart
% to snapOne) intersects with the shell wall after gauge rotation. If an
% intersection is found then the gauge positioning is not feasible and a
% tolerance is not measured. This additional check has been coded for an
% outside shell gauge placement. For nominally circular cross sections with
% small imperfections, this check is unlikely to be needed. This check does
% NOT identify any infeasible gauge placements for the sample WTST dataset.

clear;
clc;
close all

%% Nominal geometry
tower = 'ReleaseTower'; % Tower to be processed
splitPath = split(pwd,filesep);
load(fullfile(strjoin(splitPath(1:end-1),filesep), strcat('Input_', tower),strcat('S6_Mesh_', tower, '.mat')))

% Plotting used data
figure;
surf(THET_MESH, Z_MESH, R_MESH)
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
title('Cylindrical coordinates WTST', 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
grid on

%% Extracting joint locations
[~, STRAKES] = find(N2_MESH); % Finding strakes which have a meshed surface
N2_MESH = N2_MESH(STRAKES); % Extracting number of elements in available strakes
N2Cum = cumsum(N2_MESH); % Performing cumulative sum to extract values from overall meshed surface

%% Extracting nominal shell
% req = @(z, S) ((R0_TOPS(S) - R0_BOTS(S))/(Z_TOPS(S) - Z_BOTS(S))*(z - Z_BOTS(S)) + R0_BOTS(S)); % Equation defining nominal radius of shell with nominal z coordinates of joints
req2 = @(z, S, zbot, ztop) ((R0_TOPS(S) - R0_BOTS(S))/(ztop - zbot)*(z - zbot) + R0_BOTS(S)); % Equation defining nominal radius of shell with identified z coordinates of joints (more accurate)
z = Z_MESH(:,1); % Z coordinates
theta = THET_MESH(1,:); % THET coordinates
rnom = zeros(size(R_MESH,1),1); % Nominal radial coordinates
tnom = rnom; % Nominal thickness
startInd = 1;
for j = 1:length(STRAKES)
    S = STRAKES(j);
    inds = startInd:N2Cum(j); % Indicies of a particular strake
    % rnom(inds) = req(z(inds), S); % Less accurate
    rnom(inds) = req2(z(inds), S, z(inds(1)), z(inds(end)));
    tnom(inds) = THICKS(S)*ones(size(inds));
    startInd = N2Cum(j)+1;
end

rnom = rnom - tnom/2; % Correcting to obtain nominal shell midsurface

%% Entire surface evaluation
U0x = zeros(size(X_MESH));
plotting = false; % Whether plots are generated or not. Turn off for looped runs!!!

for k = 1:size(X_MESH,1)
    fprintf('Examining k = %g\n', k)
    xv = X_MESH(k,:); yv = Y_MESH(k,:);

    % Best fitting circle to ensure circular cross section has the most
    % accurate centre
    [xv, yv, thetv] = bestFit(xv, yv);
    xv = [xv xv]; yv = [yv yv]; thetv = [thetv thetv+2*pi];
    for j = 1:size(X_MESH,2)

        % Tolerance measurement for a given start point as described in
        % pseudocode
        [lgx, lgxtheta, x1, y1, xcent, ycent, perVec, snapOneInd, searchLimit] = defineCircGauge(xv, yv, thetv, rnom, tnom, j, k);
        while true
            [snapTwoPos, interPairs, intersects, snapOneX, snapOneY, gaugePointX, gaugePointY, gaugeEndX, gaugeEndY, snapTwoPosX, snapTwoPosY, searchRange] = computeShellGaugeCorrespondences(xv, yv, x1, y1, xcent, ycent, lgxtheta, rnom, searchLimit, snapOneInd, k, plotting);
            [x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, tolMeasureRange, gaugeSnapInd, snapOneInd, snapOneIndOld, endReached] = applyGaugeRotations(gaugePointX, gaugePointY, snapTwoPosX, snapTwoPosY, snapOneInd, snapOneX, snapOneY, snapTwoPos, intersects, x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, xv, yv, interPairs, plotting);
            [measureTol] = selfIntersectionCheck(x1, y1, xcent, ycent, snapOneX, snapOneY, xv, yv, rnom, j, k, snapOneIndOld, lgxtheta, plotting);
            if measureTol
                [U0Tol] = measureTolerance(x1, y1, xcent, ycent, snapOneX, snapOneY, gaugeEndX, gaugeEndY, gaugeSnapInd, tolMeasureRange, xv, yv, interPairs, rnom, k, lgx, lgxtheta, snapTwoPos, snapOneIndOld, plotting);
                tolMeasureRange(tolMeasureRange>size(X_MESH,2)) = tolMeasureRange(tolMeasureRange>size(X_MESH,2)) - size(X_MESH,2); % Allowing tolerance measurements to wrap around to start of cross section
                U0x(k, tolMeasureRange) = max([U0x(k, tolMeasureRange)' abs(U0Tol)],[],2);
            end

            if endReached
                break
            end
        end

    end
end

%% Plotting tolerance variation of U0x 
figure;
hold on
surf(THET_MESH, Z_MESH, R_MESH, U0x)
colormap("summer")
shading interp
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
hL = ylabel(cb,'U$_{0x}$', 'interpreter', 'latex'); set(hL,'Rotation',0);
hL.FontSize = 11;
set(gca,'TickLabelInterpreter','latex')
xlim([0 2*pi])
ylim([min(min(Z_MESH)) max(max(Z_MESH))])
grid on
title('U$_{0x}$', 'interpreter', 'latex')


%% FUNCTIONS
function [lgx, lgxtheta, x1, y1, xcent, ycent, perVec, snapOneInd, searchLimit] = defineCircGauge(xv, yv, thetv, rnom, tnom, j, k)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Extract coordinates of gaugeStart from shellPoints
% Initialise snapOneInd, to j representing the start of the gauge index. 
% Compute the gauge length lgx or lgθ using the nominal properties.
% Compute angle swept by the gauge θl.
% Compute normalised vector from gaugeStart to nominal cross-section centre (0,0)
% Compute coordinates of gaugeCentre assuming it lies on vector from gaugeStart to (0,0)
% Compute searchLimit, the index of the point π radians away from snapOneInd
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

r = rnom(k); t = tnom(k);

x1 = xv(j); y1 = yv(j); % Start of gauge
snapOneInd = j;

lgx = 4*sqrt(r*t); % Length of curved rod
lgxtheta = lgx/r; % Angle swept by rod

% Centroid of circle
perVec = [x1; y1]; % Perpendicular from (x1, y1)
perVec = perVec/norm(perVec); % Normalised perpendicular from (x1, y1)

xcent = x1 - r*perVec(1); % X Location of centre of arc
ycent = y1 - r*perVec(2); % X Location of centre of arc

% Not all points along the circumference are considered for snapping. Using
% points further than pi away from snapOne in an anticlockwise direction
% can cause the gauge to potentially rotate the wrong way and snap onto the
% point immediately preceeding the start point (depending on the magnitude
% of imperfections) which is incorrect. Thus points more than pi away from
% snapOne are not to be considered in the algorithm. It is only possible to
% reach a point pi away on a perfect shell. Thus the gauge is limited to an
% arc smaller than a semi circle with a radius smaller than the shell
% radius. For gauges of length lgx = 4*sqrt(rt), this is satisfied for r/t
% > (16/pi^2) i.e. r/t > 1.62 which is very likely the case for in-service
% shells. The other gauge length lg0 is limited to be smaller than r so
% will also be smaller than a semi circle.

searchLimit = find(thetv > thetv(j) & thetv < thetv(j) + pi, 1, 'last'); % < and not <= is used to allow the segment which may cross over into >pi range to be considered.

end

function [snapTwoPos, interPairsMatch, intersects, snapOneX, snapOneY, gaugePointX, gaugePointY, gaugeEndX, gaugeEndY, snapTwoPosX, snapTwoPosY, searchRange] = computeShellGaugeCorrespondences(xv, yv, x1, y1, xcent, ycent, lgxtheta, rnom, searchLimit, snapOneInd, k, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise searchRange, the indexes of shellPoints to check from snapOneInd until the end 
% of the shellPoints.
% Initialise snapOne, the coordinates of the point about which the gauge will rotate.
% Compute snapOne2EndLength, the distance between snapOne and gaugeEnd.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Finding snapOne coordinates
snapOneX = xv(snapOneInd); snapOneY = yv(snapOneInd);
r = rnom(k);

%%% CORRESPONDENCES ALONG LINE SEGMENTS ON WALL %%%
searchRange = snapOneInd:searchLimit; % Using this search range for end intersection 

% Identifying line segments where the end of the gauge can land
[gaugeEndX, gaugeEndY] = gaugePoints(x1, y1, xcent, ycent, lgxtheta, r); % End of gauge coordinates
snapOne2EndLength = sqrt((snapOneX-gaugeEndX)^2 + (snapOneY-gaugeEndY)^2); % Distance between end of gauge and SnapOne

adjacentPtLength = sqrt((gaugeEndX - snapOneX)^2 + (gaugeEndY - snapOneY)^2);
dist2AllPts = sqrt((snapOneX-xv(searchRange)).^2 +(snapOneY-yv(searchRange)).^2); % Distance to all points on the circumference from SnapOne

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find interPairs, the line segments along the shell where the end of the gauge could 
% possibly land. These are stored as pairs of indexes.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing line segments where one end is inside circle of radius
% snapOne2EndLength and origin snapOne and the other end is outside i.e.
% single intersection point.

mlPair = find(dist2AllPts(1:end-1) >= snapOne2EndLength & dist2AllPts(2:end) <= snapOne2EndLength)'; % More than r, then less than r
lmPair = find(dist2AllPts(1:end-1) <= snapOne2EndLength & dist2AllPts(2:end) >= snapOne2EndLength)'; % Less than r, then more than r
interPairs = [mlPair mlPair+1; lmPair lmPair+1;]+snapOneInd-1; % Line segments

% -------------------------------------------------------------------- %
%%% Trying to catch edge case where the ends of the line segment are both
%%% outside the above mentioned circle i.e. further than snapOne2EndLength
%%% away from snapOne but there is still an intersection (1 tangential
%%% intersection or two intersections between line and circle.)

% Possible line segments
interPairExtra = [searchRange(1:end-1); searchRange(2:end)];

% Removing line segments where both ends are inside the mentioned circle
sel = dist2AllPts(1:end-1) < snapOne2EndLength & dist2AllPts(2:end) < snapOne2EndLength;
interPairExtra(:,sel) = [];

% Start and end coordinates of remaining line segments
xvStart = xv(interPairExtra(1,:)); yvStart = yv(interPairExtra(1,:));
xvStartp1 = xv(interPairExtra(2,:)); yvStartp1 = yv(interPairExtra(2,:));

% Normalised vector along line segment
lineVector = [xvStartp1 - xvStart;yvStartp1 - yvStart;];
lineVectorLen = vecnorm([xvStartp1 - xvStart;yvStartp1 - yvStart;],2,1);
lineVector = lineVector./lineVectorLen;
a = lineVector(1,:);
b = lineVector(2,:);

% Computing coordinate of closest point on line segment (dot product of
% vector from centre of circle to line segment and vector along line
% segment must be zero)
lam = -1*(a.*(xvStart-snapOneX) + b.*(yvStart-snapOneY))./(a.^2+b.^2);
closeVect = [(xvStart-snapOneX)+a.*lam; (yvStart-snapOneY)+b.*lam];

% Computing distance to closest point on line segment from centre of circle
closeDist = sqrt(closeVect(1,:).^2 + closeVect(2,:).^2);

% Finding points where the closest point to the centre of the circle is
% also on the line segment. Removing if it is not on line segment.
sel = lam >0 & lam <= lineVectorLen;
interPairExtra = interPairExtra(:,sel)';
closeDist = closeDist(sel);

% Finding points where closest point on the line segment is closer than
% snapOne2EndLength i.e. there is an intersection possible
sel = closeDist > 0 & closeDist < snapOne2EndLength;
interPairExtra = interPairExtra(sel,:);
% -------------------------------------------------------------------- %

% Combining both solutions and removing any doubles if any
interPairs = [interPairs; interPairExtra];
interPairs = unique(interPairs, 'rows');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for pair in interPairs do
%     Solve analytically for the intersection between the line segment and a circle of radius 
%     snapOne2EndLength and origin snapOne and store in array intersects.
% end for  
% Append intersects to array snapTwoPosPoints, and gaugeEnd to gaugePoints for each 
% intersection of intersects found if any.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Solving for intersections of gauge with line segments
intersects = [];
interPairsMatch = []; % Ensuring that the interPairs have a one to one correspondence between the intersects computed as there may be more than one for each line segment
for pairInd = 1:size(interPairs,1)
    [posSols] = lineCircleIntersect(interPairs(pairInd,1), snapOneX, snapOneY, adjacentPtLength, xv, yv);
    intersects(end+1:end+size(posSols,1),:) = posSols;
    if any(posSols)
        interPairsMatch(end+1:end+size(posSols,1),:) = interPairs(pairInd,:).*ones(size(posSols));
    end
end


%%% CORRESPONDENCES WITH VERTICIES ON WALL %%%

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find snapTwoPos, the indexes of points of shellPoints with an index greater than or equal 
% to snapOneInd within snapOne2EndLength of snapOne where the gauge could rotate to.
% Extract coordinates of points along shell wall, snapTwoPosPoints, with index snapTwoPos.
% Compute points, gaugePoints, along the gauge with the same distance from snapOne as 
% snapTwoPosPoints.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Finding relevant points
searchRange = searchRange(2:end);
dist2AllPts = dist2AllPts(2:end);
distanceCheckRange = find(dist2AllPts < snapOne2EndLength);
snapTwoPos = searchRange(distanceCheckRange);

% snapTwo coordinates
snapTwoPosX = xv(snapTwoPos); snapTwoPosY = yv(snapTwoPos);

% Distance and angle enclosed from snapOne to snapTwo
adjacentPtLength = sqrt((snapTwoPosX - snapOneX).^2 + (snapTwoPosY - snapOneY).^2);
beta = 2*asin(adjacentPtLength/2/r);

% Defining points along gauge
[gaugePointX, gaugePointY] = gaugePoints(x1, y1, xcent, ycent, beta+atan2(snapOneY-ycent, snapOneX-xcent)-atan2(y1-ycent, x1-xcent), r);


% Adding end of gauge correspondences
gaugePointX = [gaugePointX gaugeEndX*ones(1,size(intersects,1))]; gaugePointY = [gaugePointY gaugeEndY*ones(1,size(intersects,1))];
snapTwoPosX = [snapTwoPosX intersects(:,1)']; snapTwoPosY = [snapTwoPosY intersects(:,2)'];


if plotting
    figure('Units', 'pixels','Position', [1920+680 458 900 500])
    tiledlayout(1,4)
    nexttile
    hold on
    plot(xv, yv, '-xk')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    plot(snapTwoPosX, snapTwoPosY, 'bo')
    % plot(xshell, yshell, 'k','Marker','x', 'MarkerFaceColor','r')
    plot(gaugePointX, gaugePointY, 'mx')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Considered snapping locations', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

end

function [x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, tolMeasureRange, gaugeSnapInd, snapOneInd, snapOneIndOld, endReached] = applyGaugeRotations(gaugePointX, gaugePointY, snapTwoPosX, snapTwoPosY, snapOneInd, snapOneX, snapOneY, snapTwoPos, intersects, x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, xv, yv, interPairs, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the angle enclosed by vectors from snapOne to gaugePoints and snapOne to 
% snapTwoPosPoints. 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing angle using cross product
av = [gaugePointX-snapOneX; gaugePointY-snapOneY; zeros(size(gaugePointX))]; bv = [snapTwoPosX-snapOneX; snapTwoPosY-snapOneY; zeros(size(gaugePointX))]; 
snapRotTwoAng = asin(cross(av, bv)./vecnorm(av,2,1)./vecnorm(bv,2,1)); % Not using absolute value of cross product to ensure the sign is correct
snapRotTwoAng = real(snapRotTwoAng(3,:)); % You are supposed to take abs, but using the final value preserves the rotation sign

% Performing obtuse angle correction accounting for angle sign
snapRotTwoAngSign = sign(snapRotTwoAng); % Storing sign of values
sel = dot(av, bv) < 0; % Where dot product is negative, angle is obtuse
snapRotTwoAng = abs(snapRotTwoAng); % Angle magnitude
snapRotTwoAng(sel) = pi - snapRotTwoAng(sel); % Obtuse correction
snapRotTwoAng = snapRotTwoAng.*snapRotTwoAngSign; % Restoring sign

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine the appropriate angle to use based on its sign and whether the inside or outside 
% of the shell is being assessed.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% In the extremely unlikely case there are two points with the exact same
% rotation angle, the one with the smaller index would be picked.

% Algorithm is to be only used with gauge placed on outside, but option of
% using inside gauge placement left available.

% side = 'inside';
side = 'outside';
if strcmp(side, 'inside')
    [~, gaugeSnapInd] = max(snapRotTwoAng);
    snapRotTwoAng = snapRotTwoAng(gaugeSnapInd);
elseif strcmp(side, 'outside')
    [~, gaugeSnapInd] = min(snapRotTwoAng);
    snapRotTwoAng = snapRotTwoAng(gaugeSnapInd);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Rotate gaugeStart and gaugeEnd about snapOne.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Rotation matrix to facilitate second snap
Rotm = [cos(snapRotTwoAng) -sin(snapRotTwoAng); sin(snapRotTwoAng) cos(snapRotTwoAng)];

% Rotating rod centre
rotated = Rotm*[xcent-snapOneX; ycent-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
xcent = rotated(1,:); ycent = rotated(2,:); % Extracting rotated coordinates

% Rotating rod start
rotated = Rotm*[x1-snapOneX; y1-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
x1 = rotated(1,:); y1 = rotated(2,:); % Extracting rotated coordinates

% Rotating gauge point - Just for plotting
rotated = Rotm*[gaugePointX-snapOneX; gaugePointY-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
gaugePointX = rotated(1,:); gaugePointY = rotated(2,:); % Extracting rotated coordinates

% Rotating gauge end - Just for plotting
rotated = Rotm*[gaugeEndX-snapOneX; gaugeEndY-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
gaugeEndX = rotated(1,:); gaugeEndY = rotated(2,:); % Extracting rotated coordinates

if plotting
    nexttile
    hold on
    plot(xv, yv, '-kx')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    plot(gaugePointX, gaugePointY, 'mx')
    plot(gaugePointX(gaugeSnapInd), gaugePointY(gaugeSnapInd), 'bo')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Snap', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine snapTwo, the second point of contact of the gauge with the shell wall.
% Determine tolMeasureRange, the indexes of points where tolerance measurements are 
% feasible.
% Determine the boolean endReached whether the end of the shell or gauge has been reached.
% if not endReached then
%     Update snapOneInd with index of snapTwo.
% end if
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Next rotation
if gaugeSnapInd <= length(snapTwoPos) % End of gauge not reached
    tolMeasureRange = snapOneInd+1:(snapOneInd + gaugeSnapInd)-1;
    snapOneIndOld = snapOneInd;
    snapOneInd = snapTwoPos(gaugeSnapInd);
    endReached = false;

% End of gauge reached
else
    tolMeasureRange = snapOneInd+1:interPairs(gaugeSnapInd-length(snapTwoPos),1);
    snapOneIndOld = snapOneInd;
    endReached = true;
end

end

function [measureTol] = selfIntersectionCheck(x1, y1, xcent, ycent, snapOneX, snapOneY, xv, yv, rnom, j, k, snapOneIndOld, lgxtheta, plotting)
numCheckPts = 100;
r = rnom(k);

% Compupting test gauge
snapOneTheta = atan2(snapOneY-ycent, snapOneX-xcent); % Position of snapOne relative to gaugeCentre
startTheta = atan2(y1-ycent, x1-xcent); % Position of start of gauge relaltive to gaugeCentre

snapOneThetaRelStart = snapOneTheta-startTheta; % Angular difference between snapOne and gaugeStart
if snapOneThetaRelStart < 0; snapOneThetaRelStart = snapOneThetaRelStart+2*pi; end % Correcting so numbers positive

% Region from gaugeStart to original snapOne
if j == snapOneIndOld % Making sure original snapOne is not same as gaugeStart
    firstRegionCheckPointsTheta = [];
else
    firstRegionCheckPointsTheta = linspace(0,snapOneThetaRelStart,numCheckPts); % Points between gaugeStart and snapOne
end

secondRegionCheckPointsTheta = linspace(snapOneThetaRelStart,lgxtheta,numCheckPts); % Points between snapOne and gaugeEnd
checkPointsTheta = [firstRegionCheckPointsTheta(2:end-1) secondRegionCheckPointsTheta(2:end-1)]; % Removing first and last values of each region to avoid end points on shell wall being considered
[gaugeCheckX, gaugeCheckY] = gaugePoints(x1, y1, xcent, ycent, checkPointsTheta, r); % Cartesian coordinates of points to consider

% Checking for intersections
[in,on] = inpolygon(gaugeCheckX,gaugeCheckY,xv,yv);

sel = ~in|on;
if any(find(~sel))
    % Intersection found
    measureTol = false;
else
    % No intersection found
    measureTol = true;
end

if plotting
    nexttile;
    hold on
    plot(xv, yv, '-kx')
    % plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    % plot(xc, yc, 'g','Marker','x', 'MarkerFaceColor','r')
    plot(gaugeCheckX(~sel), gaugeCheckY(~sel), 'ro')
    plot(gaugeCheckX(sel), gaugeCheckY(sel), 'go')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Intersection check', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

function [U0Tol] = measureTolerance(x1, y1, xcent, ycent, snapOneX, snapOneY, gaugeEndX, gaugeEndY, gaugeSnapInd, tolMeasureRange, xv, yv, interPairs, rnom, k, lgx, lgxtheta, snapTwoPos, snapOneIndOld, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute normalised vectors n, from gaugeCentre to shellPoints.
% Compute possible points on the gauge which are rnom from gaugeCentre along vectors n 
% in both opposing directions.
% Determine whether the computed points truly lie on the gauge using gaugeStart and θl and 
% remove incorrect points.
% Determine the distance between shellPoints and points on the gauge.
% Normalise distances by the gaugeLength to yield the toleranceValues.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

r = rnom(k);
% Identifying feasible tolerance measurement locations (must be within
% lgxtheta).
% For all points identified in range of snapOne and snapTwo, find points on
% the gauge which have a perpendicular line between the gauge and the shell

%%%%%%%%%%%%%%%%%%%% GAUGE FORMULATION CODE %%%%%%%%%%%%%%%%%%%%
% Normalised vector from shell to centroid
perVec = [xv(tolMeasureRange)-xcent; yv(tolMeasureRange)-ycent];
perVec = perVec./sqrt(perVec(1,:).^2+perVec(2,:).^2);

% Coordinates along arc of gauge on both halves of circle where the
% perpendicular from the gauge could originate from
xc = [xcent + r*perVec(1,:); xcent - r*perVec(1,:)];
yc = [ycent + r*perVec(2,:); ycent - r*perVec(2,:)];
tolMeasureRange = [tolMeasureRange; tolMeasureRange];

% Computing angle swept from start of gauge to points on gauge.
gaugeAngle = atan2(yc-ycent, xc-xcent); % Angular position of points on gauge relative to gaugeCentre
gaugeStartAngle = atan2(y1-ycent, x1-xcent); % Angular position of gaugeStart relative to gaugeCentre
arcAngle = mod(gaugeAngle-gaugeStartAngle,2*pi); % Clockwise angle swept from gaugeStart to points on gauge

% Identifying feasible coordinates of points on gauge based on whether the
% perpendicular lands within the range of the gauge.
sel = arcAngle <= lgxtheta;
xc = xc(sel); yc = yc(sel);
tolMeasureRange = tolMeasureRange(sel);

% %%%%%%%%%%%%%%%%%%%% GAUGE FORMULATION CODE %%%%%%%%%%%%%%%%%%%%

% Extracting dimple tolerance - There may be points which satisfy all
% conditions but do not have a direct line of sight from the gauge to the
% point on the shell wall due to the shell wall itself obstructing the line
% of sight. Tolerance measures are still taken in these cases.

radialDist = sqrt((xc-xv(tolMeasureRange)').^2 + (yc-yv(tolMeasureRange)').^2);
U0Tol = radialDist/lgx;

% Plotting
if plotting
    nexttile
    hold on
    plot(xv, yv, '-kx')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xc, yc, 'ok')
    plot([xv(tolMeasureRange); xc'], [yv(tolMeasureRange); yc'], 'g')

    % Plotting appropriate points along gauge
    if gaugeSnapInd > length(snapTwoPos)
        alongShellIndsPlotting = snapOneIndOld:interPairs(gaugeSnapInd-length(snapTwoPos),1);
    else
        alongShellIndsPlotting = snapOneIndOld:(snapOneIndOld + gaugeSnapInd);
    end

    perVec = [[xv(snapOneIndOld) xv(alongShellIndsPlotting) gaugeEndX]-xcent; [yv(snapOneIndOld) yv(alongShellIndsPlotting) gaugeEndY]-ycent]; % Perpendicular to rod vector
    perVec = perVec./sqrt(perVec(1,:).^2+perVec(2,:).^2); % Normalised Perpendicular to rod vector

    % Coordinates along arc
    gaugeCheckX = [xcent + r*perVec(1,:); xcent - r*perVec(1,:)];
    gaugeCheckY = [ycent + r*perVec(2,:); ycent - r*perVec(2,:)];

    gaugeAngle = atan2(gaugeCheckY-ycent, gaugeCheckX-xcent); % Angular position of points on gauge relative to gaugeCentre
    gaugeStartAngle = atan2(y1-ycent, x1-xcent); % Angular position of gaugeStart relative to gaugeCentre
    arcAngle = mod(gaugeAngle-gaugeStartAngle,2*pi); % Clockwise angle swept from gaugeStart to points on gauge

    % Identifying places where the perpendiculars do not land on the gauge
    sel = arcAngle <= lgxtheta;
    gaugeCheckX = gaugeCheckX(sel); gaugeCheckY = gaugeCheckY(sel);

    plot(gaugeCheckX, gaugeCheckY, 'bx')
    plot(xcent, ycent, 'ko')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Tolerance extraction', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

end

function [gaugePointX, gaugePointY] = gaugePoints(x1, y1, xcent, ycent, beta, r)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the angular position corresponding to beta.
% Convert the angular position corresponding to beta into Cartesian coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Defining points along gauge
startAng = atan2(y1-ycent, x1-xcent); % Angular position of (x1, y1) relative to gaugeCentre
endAng = startAng + beta; % Angular position of gaugePoints relative to gaugeCentre

% Coordinates of gaugePoints
gaugePointX = r*cos(endAng) + xcent;
gaugePointY = r*sin(endAng) + ycent;
end

% Solving for the intersection of shell wall and end of gauge
function [sols] = lineCircleIntersect(startInd, snapOneX, snapOneY, adjacentPtLength, xv, yv)

roundingTol = 1e-9;

% Defining circle with centre (snapOneX, snapOneY) and radius
% adjacentPtLength which relates to where the end of the gauge could lie.
s1x = snapOneX;
s1y = snapOneY;
a = adjacentPtLength;

% Defining points on line segment where circle intersects
xvStart = xv(startInd); yvStart = yv(startInd);
xvStartp1 = xv(startInd+1); yvStartp1 = yv(startInd+1);

% Computing angle of inclination of line segment - used to avoid edge cases
% such as infinite gradient of line segment
wallAng = atan2(yvStartp1 - yvStart, xvStartp1 - xvStart);
if wallAng < 0
    wallAng = wallAng + 2*pi;
end

% Computing rotational angle to rotate line onto a 45 degree diagonal
rotAng = wallAng - pi/4;

% Rotation matrix to rotate line onto a 45 degree diagonal
rotAng = -rotAng;
Rotm = [cos(rotAng) -sin(rotAng); sin(rotAng) cos(rotAng)];

% Rotating key points
rotated = Rotm*[xvStartp1-xvStart s1x-xvStart; yvStartp1-yvStart s1y-yvStart]+[xvStart; yvStart]; % Performing rotation
xvStartp1 = rotated(1,1); yvStartp1 = rotated(2,1);
s1x = rotated(1,2); s1y = rotated(2,2);

% Line defining adjacent points on shell
m = (yvStartp1 - yvStart)/(xvStartp1 - xvStart);
c = yvStartp1-m*xvStartp1;

% Solving for the intersection between line and circle
% syms m x c s1x s1y a
% eqn = s1x + sqrt(a^2 - (m*x+c - s1y)^2) - x == 0;
% S = solve(eqn,x,"ReturnConditions",true);
xEnd1 = (s1x - c*m + m*s1y + (a^2*m^2 + a^2 - c^2 - 2*c*m*s1x + 2*c*s1y - m^2*s1x^2 + 2*m*s1x*s1y - s1y^2)^(1/2))/(m^2 + 1);
xEnd2 = (s1x - c*m + m*s1y - (a^2*m^2 + a^2 - c^2 - 2*c*m*s1x + 2*c*s1y - m^2*s1x^2 + 2*m*s1x*s1y - s1y^2)^(1/2))/(m^2 + 1);
yEnd1 = m*xEnd1+c;
yEnd2 = m*xEnd2+c;

% Packaging solutions
sols = [xEnd1 yEnd1; xEnd2 yEnd2];

% Checking if solutions are in the domain of line segment (as there will always be two solutions)
sol1Feas = xEnd1 >= min([xvStart xvStartp1])-roundingTol & xEnd1 <= max([xvStart xvStartp1])+roundingTol & yEnd1 >= min([yvStart yvStartp1])-roundingTol & yEnd1 <= max([yvStart yvStartp1])+roundingTol;
sol2Feas = xEnd2 >= min([xvStart xvStartp1])-roundingTol & xEnd2 <= max([xvStart xvStartp1])+roundingTol & yEnd2 >= min([yvStart yvStartp1])-roundingTol & yEnd2 <= max([yvStart yvStartp1])+roundingTol;

% Keeping feasible solutions (ussually only one will remain)
sols = sols([sol1Feas sol2Feas],:);

% Rotation matrix to rotate solution back onto original coordinate system
rotAng = -rotAng;
Rotm = [cos(rotAng) -sin(rotAng); sin(rotAng) cos(rotAng)];

% Rotating key points
rotated = Rotm*(sols'-[xvStart;yvStart])+[xvStart; yvStart]; % Performing rotation
sols = rotated';

% Removing imaginary solutions if any (there should not be any)
if ~isreal(sols)
    sols = [];
end

end

% Best fit circle to data
function [xv, yv, thetv] = bestFit(xv, yv)
xm = [0 0]; % Initial guess

options = optimoptions('lsqnonlin');
options.MaxFunctionEvaluations = 500000;
options.FunctionTolerance = 1e-12;
options.StepTolerance = 1e-12;
options.MaxIterations = 5000;

f = @(b) circlefun(b, xv', yv');
xs = lsqnonlin(f,xm,[],[], [], [],[],[],[],options);

% Reiorientating geometry
xv = xv-xs(1);
yv = yv-xs(2);

% Compputing circumferential coordinate
thetv = atan2(yv, xv);

% Ensuring theta starts at zero contiuously and increases
thetv = thetv - thetv(1); 
thetv(thetv < 0) = thetv(thetv < 0) + 2*pi;
end

% Evaluation metric for circle fitting
function [val] = circlefun(b, X, Y)

coords = [X-b(1) Y-b(2)]; % Shifting origin
val = sqrt(coords(:,1).^2 + coords(:,2).^2); % Quantity to minimise

end