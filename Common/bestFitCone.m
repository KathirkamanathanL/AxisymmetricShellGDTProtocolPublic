function [D] = bestFitCone(X, Y, Z, minZ, tower, nRandInds)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 21.57 on 28/04/2025


% Function returns best fitting cone parameters.

% nRandInds is the number of points to consider for the best fitting as the
% entire dataset is not neccessary.

%% Loading nominal properties
run(fullfile(strcat('Input_', tower),strcat('S0_Data_', tower)))

%% Downsampling data
X = reshape(X, 1, length(X));
Y = reshape(Y, 1, length(X));
Z = reshape(Z, 1, length(X));

randInds = randsample(length(X),nRandInds); % Random idexes

Xrand = X(randInds); Yrand = Y(randInds); Zrand = Z(randInds);
xm = [Xrand; Yrand; Zrand]';

%% Best Fitting - Cone
options = optimoptions('lsqnonlin');
options.MaxFunctionEvaluations = 500000;
options.FunctionTolerance = 1e-14;
options.StepTolerance = 1e-14;
options.MaxIterations = 5000;

% Initial guess [X1 Y1 rotAng1 rotAng2 BotRadius TopRadius]
b0 = [0 0 0 0 R0_BOTS(1) R0_BOTS(1)]/1e3;

% Shifting bottom of tower to 0 m elevation
Z = Z - minZ;

% Solving for best fit parameters
tic
f = @(b) coneFit(b, xm);
D = lsqnonlin(f,b0,[],[], [], [],[],[],[],options);
toc

% Shifting tower back up
Z = Z + minZ;

end

%%
% Cone function similar to the formulation in de Vries PhD
function [perpImp] = coneFit(b, x)

% b(1) = X1
% b(2) = Y1
% b(3) = rotAng1
% b(4) = rotAng2
% b(5) = Botradius
% b(6) = Topradius

X = x(:,1); % Raw X coords
Y = x(:,2); % Raw Y coords
Z = x(:,3); % Raw Z coords
H = max(Z);

coords = [X-b(1) Y-b(2) Z]; % Shifting origin
Rx = [1 0 0; 0 cos(b(3)) -sin(b(3)); 0 sin(b(3)) cos(b(3))]; % Rotation about x axis
Ry = [cos(b(4)) 0 sin(b(4)); 0 1 0; -sin(b(4)) 0 cos(b(4))]; % Rotation about Y axis

coords = coords*Rx*Ry; % Rotating coordinates

R = (b(5) - b(6))/(H -0)*(H - Z)+b(6);
beta = atan2((b(5) - b(6)), H);

radialImp = R - sqrt((coords(:,1)).^2 + (coords(:,2)).^2);
perpImp = abs(radialImp.*cos(beta));

end