% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 22.36 on 06/05/2025


% Code measures the unintended eccentricity tolerance

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
R_MESH = R_MESH - tnom/2; % Correcting to obtain real shell midsurface

%% Computing Ue
eaM = zeros(length(STRAKES)-1, length(theta)); % Array to store unintended eccentricity [mm]
etotM = zeros(length(STRAKES)-1, length(theta)); % Array to store total eccentricity [mm]
UeM = zeros(length(STRAKES)-1, length(theta)); % Array to store Unintended eccentricity tolerance [-]

UeHeatMap = zeros(size(Z_MESH)); % Unintended eccentricity tolerance heat map [-]

for j = 1:length(STRAKES)-1
    indBot = N2Cum(j)+1; % Bottom of strake above joint ind
    indTop = N2Cum(j); % Top of strake below joint ind

    eint = rnom(indTop) - rnom(indBot); % Intended eccentricity
    etot = R_MESH(indTop,:) - R_MESH(indBot,:); % Total eccentricity
    ea = etot-eint; % Unintended eccentricity

    eaM(j,1:length(theta)) = ea; % Storing unintended eccentricty
    etotM(j,1:length(theta)) = etot; % Storing Total eccentricty

    tav = mean([tnom(indBot) tnom(indTop)]); % Average thickness at a joint

    Ue = ea/tav; % Ue tolerance
    UeM(j,1:length(theta)) = Ue; % Storing tolerance
    
    UeHeatMap(indTop:indBot,1:length(theta)) = [Ue; Ue]; % Storing data for heat map
end

%% Plotting unintended eccentricy tolerance [-]
Uemrs = reshape(UeM,size(UeM,1)*size(UeM,2),1);

LW = 1.1; % Line width
classA = 0.14;
classB = 0.2;
classC = 0.3;
varName = 'U_e';

figure;
hold on
h = histogram(abs(Uemrs),'normalization','probability');
p1 = xline(classA,'color',[0.8500 0.3250 0.0980],'Linestyle','--', 'Linewidth',LW);
p2 = xline(classB,'color',[0.9290 0.6940 0.1250],'Linestyle','--', 'Linewidth',LW);
p3 = xline(classC,'color',[0.4940 0.1840 0.5560],'Linestyle','--', 'Linewidth',LW);
xlabel(sprintf('$|%s|$ [-]',varName),'interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title(sprintf('Histogram of unintended eccentricity tolerance $|%s|$', varName),'interpreter','latex')
legend([p1, p2, p3],{'FTQ Class A', 'FTQ Class B', 'FTQ Class C'},'interpreter','latex')
xlim([0 max(max(h.BinEdges),classC)*1.1])
set(gca,'TickLabelInterpreter','latex')
grid on

figure;
hold on
h = cdfplot(abs(Uemrs));
h.LineWidth = LW;
p1 = xline(classA,'color',[0.8500 0.3250 0.0980],'Linestyle','--', 'Linewidth',LW);
p2 = xline(classB,'color',[0.9290 0.6940 0.1250],'Linestyle','--', 'Linewidth',LW);
p3 = xline(classC,'color',[0.4940 0.1840 0.5560],'Linestyle','--', 'Linewidth',LW);
xlabel(sprintf('$|%s|$ [-]',varName),'interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title(sprintf('ECDF of unintended eccentricity tolerance $|%s|$', varName),'interpreter','latex')
legend([p1, p2, p3],{'FTQ Class A', 'FTQ Class B', 'FTQ Class C'},'interpreter','latex', 'Location','southeast')
xlim([0 classC*1.1])
set(gca,'TickLabelInterpreter','latex')
grid on

figure;
surf(THET_MESH, Z_MESH, abs(UeHeatMap))
shading interp
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$|$Ue$|$ [-]', 'Interpreter','latex')
xlim([0 2*pi])
ylim([min(z) max(z)])
view([-486.1125 73.2574])
set(gca,'TickLabelInterpreter','latex')
grid on
