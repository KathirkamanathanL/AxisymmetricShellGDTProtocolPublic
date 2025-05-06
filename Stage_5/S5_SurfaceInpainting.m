function [THETg, Zg, RBIG] = S5_SurfaceInpainting(THETg, Zg, RBIG, JOINTCUTOFFBOT, JOINTCUTOFFTOP)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 15.45 on 29/04/2025


% Function inpaints gaps left by joint removal

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Add top2JointCentre to the top vertical coordinate to yield zvTop at the joint centre.
% Extrapolate the radial coordinate rvTop to the top of the meridian, zvTop.
% 
% Take bot2JointCentre from the bottom vertical coordinate to yield zvBot at the joint centre.
% Extrapolate the radial coordinate rvBot to the bottom of the meridian, zvBot.
% 
% Create new arrays rv = [rvBot rv rvTop], zv = [zvBot zv zvTop].
% Correct vertical coordinates zv such that they start at zero.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

%%%%%%%%%%% TOP OF STRAKE EXTRAPOLATION %%%%%%%%%%%
zvTop = Zg(end) + JOINTCUTOFFTOP;
rvTop = RBIG(end,:) + JOINTCUTOFFTOP*(RBIG(end,:) - RBIG(end-1,:))./(Zg(end,:) - Zg(end-1,:)); % Extrapolating top of strake

%%%%%%%%%%% BOTTOM OF STRAKE EXTRAPOLATION %%%%%%%%%%%
zvBot = -JOINTCUTOFFBOT; 
rvBot = RBIG(1,:) + (-JOINTCUTOFFBOT)'*(RBIG(2,:) - RBIG(1,:))./(Zg(2,:) - Zg(1,:)); % Extrapolating bottom of strake

%%%%%%%%%%% CREATING NEW ARRAYS %%%%%%%%%%%
% Vertical coordinates
Zg = [zvBot*ones(1,size(THETg,2)); Zg; zvTop*ones(1,size(THETg,2))]; 
Zg = Zg + JOINTCUTOFFBOT;

% Radial coordinates
RBIG = [rvBot; RBIG; rvTop];

% Circumferential coordiantes
THETg = [THETg(1,:).*ones(1,size(THETg,2)); THETg; THETg(1,:).*ones(1,size(THETg,2)); ];

end
