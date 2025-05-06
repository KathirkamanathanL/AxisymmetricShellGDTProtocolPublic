function [X_MESH, Y_MESH, Z_MESH, THET_MESH, R_MESH, N1_MESH, N2_MESH] = S6_ProjectionToArbitraryMeshes(CANS, cloudStrakes, Z_BOTS, Z_TOPS, R0_BOTS, R0_TOPS, THICKS)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 20.58 on 02/05/2025


% Function generates Koiter ellipse compliant mesh

X_MESH = []; Y_MESH = []; Z_MESH = []; THET_MESH = []; R_MESH = [];

%% Computing Koiter ellipse meshing requirements of each strake
circ_mesh_n = 0;
merd_mesh_m = zeros(size(Z_BOTS));
e = 1.0;
for S = 1:length(Z_BOTS)
    % Assigning zero elements for not meshed cans
    if ismember(S, cloudStrakes)
        rb = R0_BOTS(S); rt = R0_TOPS(S); H = Z_TOPS(S) - Z_BOTS(S); t = THICKS(S);
        beta = atan(abs(rb-rt)/H);
        rhoAv = (rb + rt)/2.0/cos(beta);
        LL = H/cos(beta);
        circ_mesh_n = max(circ_mesh_n, 18.18*e*cos(beta)*sqrt(rhoAv/t)); % Minimum number of circumferential elements
        merd_mesh_m(S) = ceil((5.79*e*LL/sqrt(rhoAv*t))); % Minimum number of meridional elements

        merd_mesh_m(S) = max([merd_mesh_m(S) 2]); % Ensuring each can has at least two elements
    else
        merd_mesh_m(S) = 0;
    end
end

merd_mesh_m = merd_mesh_m + 1; % Increasing by one to get number of nodes
circ_mesh_n = ceil(circ_mesh_n); % Do not need to increase by one as in circumferential direction, num elem = num nodes

%% Looping through cans
for S = cloudStrakes

    % Setting target mesh size of the current strake
    targN = circ_mesh_n;
    targM = merd_mesh_m(S);

    % Storing the spacing of the actual mesh
    N1_MESH(S) = targN;
    N2_MESH(S) = targM;

    % Extracting relevant part of gridded surface
    X3D = CANS.(sprintf('S%g',S)).X3D;
    Y3D = CANS.(sprintf('S%g',S)).Y3D;
    Z3D = CANS.(sprintf('S%g',S)).Z3D;

    theta = CANS.(sprintf('S%g',S)).THETg(1,:);
    THETg = ones(size(X3D)).*theta;
    Rg = sqrt(X3D.^2 + Y3D.^2);
    Zg = Z3D;

    % Adding missing column back in to ensure uniform mesh
    THETg = [THETg 2*pi*ones(size(THETg,1),1)];
    Rg = [Rg Rg(:,1)];
    Zg = [Zg Zg(:,1)];

    %%%%%%%%%%%% SCATTERED INTERPOLANT %%%%%%%%%%%%
    % Creating a meshgrid with desired spacing
    theta = linspace(min(THETg(1,:)), max(THETg(1,:)), targN+1);  theta = theta(1:end-1);
    z = linspace(min(Zg(:,1)), max(Zg(:,1)), targM);
    [THETgd, Zgd] = meshgrid(theta,z);

    % Performing interpolation
    F = scatteredInterpolant(reshape(THETg, size(THETg,1)*size(THETg,2), 1), reshape(Zg, size(THETg,1)*size(THETg,2), 1), reshape(Rg, size(THETg,1)*size(THETg,2), 1));
    Rgd = F(THETgd, Zgd);
    %%%%%%%%%%%% SCATTERED INTERPOLANT %%%%%%%%%%%%

    % Converting points back into cartesian coordinates
    Xgd = Rgd.*cos(THETgd); Ygd = Rgd.*sin(THETgd);

    % Storing mesh coordinates
    X_MESH(end+1:end+size(Xgd, 1),:) = Xgd;
    Y_MESH(end+1:end+size(Xgd, 1),:) = Ygd;
    Z_MESH(end+1:end+size(Xgd, 1),:) = Zgd;
    THET_MESH(end+1:end+size(Xgd, 1),:) = THETgd;
    R_MESH(end+1:end+size(Xgd, 1),:) = Rgd;

    % Plotting mesh of a particular can
    figure;
    subplot(1,2,1)
    hold on
    surf(THETg, Zg, Rg, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none')
    surf(THETgd, Zgd, Rgd)
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    grid on
    title(sprintf('Projected mesh of strake at %g - %g m',Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gcf,'Position',[297 275 863 516])
    set(gca, 'View',[-16.7552 56.8675])

    subplot(1,2,2)
    hold on
    surf(X3D, Y3D, Z3D, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none')
    surf(Xgd, Ygd, Zgd)
    xlabel('$x$ [mm]','interpreter','latex')
    ylabel('$y$ [mm]','interpreter','latex')
    zlabel('$z$ [mm]','interpreter','latex')
    grid on
    title(sprintf('Projected mesh of strake at %g - %g m',Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gcf,'Position',[297 275 863 516])
    set(gca, 'View',[-16.7552 56.8675])
    axis equal tight

end

% Plotting mesh of entire shell
figure
hold on
surf(THET_MESH, Z_MESH, R_MESH)
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$z$ [mm]','interpreter','latex')
zlabel('$\rho$ [mm]','interpreter','latex')
grid on
title('Projected mesh of shell','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])

figure
hold on
surf(X_MESH, Y_MESH, Z_MESH)
xlabel('$x$ [mm]','interpreter','latex')
ylabel('$y$ [mm]','interpreter','latex')
zlabel('$z$ [mm]','interpreter','latex')
grid on
title('Projected mesh of shell','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[297 275 863 516])
set(gca, 'View',[-16.7552 56.8675])
axis equal tight