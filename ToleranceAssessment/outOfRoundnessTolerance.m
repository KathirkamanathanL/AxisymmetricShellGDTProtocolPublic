% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 22.10 on 06/05/2025


% Code measures the out-of-roundness tolerance

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

%% Computing Out of roundness

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for point in (xv, yv) do
%       Compute distance between point and every other point in (xv, yv).
%       Determine maximum distance and store in a 1D array D. This corresponds to the diameter 
%       measured from the point.
% end for
% Extract the smallest and largest diameter in D and store as dmin and dmax respectively.
% Compute the out-of-roundness Ur of the cross-section using dmin, dmax and dnom.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

D = zeros(size(THET_MESH)); % Creating empty array to store actual diameters

% Working out diameter when starting from every point in the shell
for j = 1:size(THET_MESH,1)
    for k = 1:size(THET_MESH,2)
        D(j,k) = max(sqrt((X_MESH(j,k) - X_MESH(j,:)).^2 + (Y_MESH(j,k) - Y_MESH(j,:)).^2));
    end
end

Ur = (max(D,[],2) - min(D,[],2))./(2*rnom); % Computing out of roundness

% Plotting histogram of Out of Roundness tolerance Ur
figure;
hold on
h = histogram(Ur,'normalization','probability');
p1 = xline(0.007,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.8500 0.3250 0.0980]);
p2 = xline(0.010,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.9290 0.6940 0.1250]);
p3 = xline(0.015,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.4940 0.1840 0.5560]);
xlabel('$U_r$ [-]','interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title('Histogram of Out-of-roundness tolerance $U_r$','interpreter','latex')
legend([p1, p2, p3],{'FTQ Class A', 'FTQ Class B', 'FTQ Class C'},'interpreter','latex')
xlim([0 max(max(h.BinEdges),0.015)*1.1])
set(gca,'TickLabelInterpreter','latex')
grid on


% Plotting variation of Out of Roundness tolerance Ur with meridional height z
figure;
hold on
plot(Ur, z, 'Linewidth',1.1)
p1 = xline(0.007,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.8500 0.3250 0.0980]);
p2 = xline(0.010,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.9290 0.6940 0.1250]);
p3 = xline(0.015,'Linestyle', '--', 'Linewidth', 1.1, 'Color', [0.4940 0.1840 0.5560]);
xlabel('$U_r$ [-]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
title('Variation of Out-of-roundness tolerance $U_r$ with height','interpreter','latex')
legend([p1, p2, p3],{'FTQ Class A', 'FTQ Class B', 'FTQ Class C'},'interpreter','latex')
xlim([0 max(max(Ur),0.015)*1.1])
set(gca,'TickLabelInterpreter','latex')
grid on