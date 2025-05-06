function [TT_REM, ZZ_REM, RR_REM, TT_OUT_BELOW, ZZ_OUT_BELOW, RR_OUT_BELOW, TT_OUT_ABOVE, ZZ_OUT_ABOVE, RR_OUT_ABOVE, DISTPQ] = S3_SystematicAndRandomOutlierRemoval(TT, ZZ, RR, RMEAN,  tbinSize, zbinSize, cutOff)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 17.52 on 02/05/2025


% Function removes systematic and random outliers

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Scale circumferential coordinate by rmean.
% Extend points circumferentially.
% Compute circumferential and vertical coordinates used to define sliding windows
% Initialise an array distpq to store the maximum projected distance of every point to its plane.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Ensuring the z coordinate starts from zero
minZZ = min(ZZ);
ZZ = ZZ - minZZ;

% Ensuring theta varies from 0 to 2pi
TT(TT < 0) = TT(TT < 0) + 2*pi;

% Scaling theta to a reasonable value
TT = TT *RMEAN;

% Extending point cloud from 0 radians to -tbinSize/RMEAN
rightsel = TT > max(TT) - tbinSize;
rightOriginalInds = find(rightsel);
rightAddedInds = length(TT)+1:length(TT)+length(rightOriginalInds);
TT(end+1:end+length(TT(rightsel))) = TT(rightsel) - 2*pi*RMEAN;
ZZ(end+1:end+length(TT(rightsel))) = ZZ(rightsel);
RR(end+1:end+length(TT(rightsel))) = RR(rightsel);

% Extending point cloud from 2*pi radians to 2*pi+tbinSize/RMEAN
leftsel = TT < tbinSize;
leftOriginalInds = find(leftsel);
leftAddedInds = length(TT)+1:length(TT)+length(leftOriginalInds);
TT(end+1:end+length(TT(leftsel))) = TT(leftsel) + 2*pi*RMEAN;
ZZ(end+1:end+length(TT(leftsel))) = ZZ(leftsel);
RR(end+1:end+length(TT(leftsel))) = RR(leftsel);

% Creating bin coordinates for z and theta with an assumed bin size of zbinSize and tbinSize
z = linspace(0, max(ZZ), ceil(max(ZZ)/zbinSize)*3);
theta = linspace(min(TT), max(TT), ceil(max(TT)/tbinSize)*3);

% Array storing distances between pca plane and point
DISTPQ = zeros(size(TT));

%% Looping through sliding windows

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for window in slindingWindows do
%     Identify points which belong to the window
%     Perform PCA to identify the surface normal
%     Select a point at random through which the plane passes.
%     Project points onto the plane and compute the distance to the plane.
%     Find the point with the median projected signed distance to the plane and redefine the plane to pass through this point.
%     Compute the new projected distance of points to the new plane.
%     Update distpq for all points in the window with their new projected distance provided it is larger than their existing value in distpq.
% end for
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

for k = 1:length(z) - 3
    for j = 1:length(theta) - 3

        % Isolating points in bin
        sel = ZZ > z(k) & ZZ < z(k+3) & TT > theta(j) & TT < theta(j+3);
        TTT = TT(sel);
        ZZZ = ZZ(sel);
        RRR = RR(sel);
        
        % Seeing if there are enough points for PCA
        if length(TTT) <= 6
            fprintf('Inadequate data in j = %g, k = %g\n', j, k)
            continue
        end

        % Performing pca and extracting the normal vector n
        n = pca([TTT' ZZZ' RRR']);
        n = n(:,3);

        % Finding vector between a random point (which is set as a point on the plane) and all other points in bin
        distpq = [TTT' ZZZ' RRR'] - [TTT(1) ZZZ(1) RRR(1)]; % First point chosen as random start point which will be corrected in a few lines
        
        % Computing projected signed distance between all the points and the plane
        distpq = distpq*n;

        % Finding the location of the point which is a median signed distance from
        % the plane. This should be a point which is roughly in the middle
        % of all points when considering the perpendicular direction to the
        % plane
        mm = find(min(abs(distpq - median(distpq))) == abs(distpq - median(distpq)),1);

        % Recomputing the distances after shifting the assumed point on the plane
        distpq = [TTT' ZZZ' RRR'] - [TTT(mm) ZZZ(mm) RRR(mm)];
        distpq = distpq*n;

        % Plotting
        % figure
        % hold on
        % scatter3(TTT(m), ZZZ(m), RRR(m), 200, '.b')
        % scatter3(TTT, ZZZ, RRR, 10, '.r')
        % scatter3(TTT((abs(distpq) > cutOff)), ZZZ((abs(distpq) > cutOff)), RRR((abs(distpq) > cutOff)), 100, '.g')

        % Storing the computed distances in DISTPQ
        compVals = DISTPQ(sel); % Values already calculated
        compValsSel = abs(compVals) > abs(distpq'); % Where old vals bigger than new vals
        distpq((compValsSel)) = compVals(compValsSel); % Replace new vals vector with larger old vals
        DISTPQ(sel) = distpq; % Update the big vector 
    end
end

%% Identifying outliers

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Restore circumferential coordinate to its value prior to scaling.
% Remove points added to data when extending circumferentially.
% Identify outlierFreePoints with a distpq less than cutoff.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

TT = TT /RMEAN; % Rescaling
ZZ = ZZ + minZZ;

% Updating values within 0 and 2 pi with values of points outside the range
sel = abs(DISTPQ(rightAddedInds)) > abs(DISTPQ(rightOriginalInds));
DISTPQ(rightOriginalInds(sel)) = DISTPQ(rightAddedInds(sel));

sel = abs(DISTPQ(leftAddedInds)) > abs(DISTPQ(leftOriginalInds));
DISTPQ(leftOriginalInds(sel)) = DISTPQ(leftAddedInds(sel));

% Removing added points
sel = TT < 0;
DISTPQ(sel) = []; TT(sel) = []; ZZ(sel) = []; RR(sel) = [];

sel = TT > 2*pi;
DISTPQ(sel) = []; TT(sel) = []; ZZ(sel) = []; RR(sel) = [];

% Outlier free points
outlierFreeSel = (abs(DISTPQ) < cutOff);
TT_REM = TT(outlierFreeSel);
ZZ_REM = ZZ(outlierFreeSel);
RR_REM = RR(outlierFreeSel);

% Outliers above cloud
outlierFreeSelAbove = DISTPQ > cutOff;
TT_OUT_ABOVE = TT(outlierFreeSelAbove);
ZZ_OUT_ABOVE = ZZ(outlierFreeSelAbove);
RR_OUT_ABOVE = RR(outlierFreeSelAbove);

% Outliers below cloud
outlierFreeSelBelow = DISTPQ < -cutOff;
TT_OUT_BELOW = TT(outlierFreeSelBelow);
ZZ_OUT_BELOW = ZZ(outlierFreeSelBelow);
RR_OUT_BELOW = RR(outlierFreeSelBelow);

% Plotting
figure
hold on
scatter3(TT_REM, ZZ_REM, RR_REM, 100, '.k')
scatter3(TT_OUT_BELOW, ZZ_OUT_BELOW, RR_OUT_BELOW, 100, '.r')
scatter3(TT_OUT_ABOVE, ZZ_OUT_ABOVE, RR_OUT_ABOVE, 100, '.b')
grid on
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])

% Reshaping output
TT_REM = reshape(TT_REM, [length(TT_REM), 1]);
ZZ_REM = reshape(ZZ_REM, [length(ZZ_REM), 1]);
RR_REM = reshape(RR_REM, [length(RR_REM), 1]);


end