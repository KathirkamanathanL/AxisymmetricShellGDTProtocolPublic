classdef S1_CLOUD

    % Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
    % of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
    % Energy.

    % Copyright under a BSD 3-Clause License, see
    % https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

    % Last modified at 16.34 on 29/04/2025


    % Class used to perform actions on individual point clouds

    properties
        ID % Identifier for cloud
        X; Y; Z % XYZ coordinates
        R; T % Cylindrical coordinates
        XNEW; YNEW; ZNEW % Registered coordinates
        normals % Point cloud normals
    end


    methods

        function obj = S1_CLOUD(fileName, ID)

            fileID = fopen(fileName);
            B = fread(fileID,'double');
            fclose(fileID);
            B = B(2:end);
            B = reshape(B,[3, length(B)/3])';

            obj.X = B(:,1);
            obj.Y = B(:,2);
            obj.Z = B(:,3);
            obj.ID = ID;

        end

        % Downsamples the cloud
        function obj = CLOUDDOWNRND(obj,nperc)

            ptCloudIn = pointCloud([obj.X obj.Y obj.Z]);
            ptCloudOut = pcdownsample(ptCloudIn,'random',nperc);
            obj.X = ptCloudOut.Location(:,1); obj.Y = ptCloudOut.Location(:,2); obj.Z = ptCloudOut.Location(:,3);

        end

        % Limits Z range of cloud
        function obj = CLOUDZLIM(obj,minZ,maxZ)

            sel = obj.Z >= minZ & obj.Z <= maxZ;
            obj.X = obj.X(sel);
            obj.Y = obj.Y(sel);
            obj.Z = obj.Z(sel);

        end

        % Returns the R and T cylindrical coordinates
        function obj = CLOUDRT(obj, X, Y)

            obj.R = sqrt(X.^2 + Y.^2);
            obj.T = atan2(Y,X);

        end

        % Returns the registered cloud
        function obj = CLOUDMERGE(obj, C)

            X0 = obj.X; Y0 = obj.Y; Z0 = obj.Z;
            nopt = 6;
            Ci = C(nopt*obj.ID-(nopt-1):nopt*obj.ID);
            Rxi = [1 0 0; 0 cos(Ci(1)) -sin(Ci(1)); 0 sin(Ci(1)) cos(Ci(1))]; % Rotation about x axis
            Ryi = [cos(Ci(2)) 0 sin(Ci(2)); 0 1 0; -sin(Ci(2)) 0 cos(Ci(2))]; % Rotation about Y axis
            Rzi = [cos(Ci(3)) -sin(Ci(3)) 0; sin(Ci(3)) cos(Ci(3)) 0; 0 0 1]; % Rotation about Y axis

            fc = Rzi*Ryi*Rxi*[X0'; Y0'; Z0'] + [Ci(4); Ci(5); Ci(6)];

            obj.XNEW = fc(1,:);
            obj.YNEW = fc(2,:);
            obj.ZNEW = fc(3,:);

        end

        % Returns the iteratively registered cloud
        function obj = CLOUDMERGEITER(obj, C, iuse)
            nopt = 6;
            X0 = obj.X; Y0 = obj.Y; Z0 = obj.Z;

            for solNum = 1:iuse
                X0 = reshape(X0, length(X0),1);
                Y0 = reshape(Y0, length(Y0),1);
                Z0 = reshape(Z0, length(Z0),1);

                Ci = C(nopt*obj.ID-(nopt-1):nopt*obj.ID, solNum);
                Rxi = [1 0 0; 0 cos(Ci(1)) -sin(Ci(1)); 0 sin(Ci(1)) cos(Ci(1))]; % Rotation about x axis
                Ryi = [cos(Ci(2)) 0 sin(Ci(2)); 0 1 0; -sin(Ci(2)) 0 cos(Ci(2))]; % Rotation about Y axis
                Rzi = [cos(Ci(3)) -sin(Ci(3)) 0; sin(Ci(3)) cos(Ci(3)) 0; 0 0 1]; % Rotation about Z axis

                fc = Rzi*Ryi*Rxi*[X0'; Y0'; Z0'] + [Ci(4); Ci(5); Ci(6)];

                X0 = fc(1,:);
                Y0 = fc(2,:);
                Z0 = fc(3,:);
            end

            obj.XNEW = fc(1,:);
            obj.YNEW = fc(2,:);
            obj.ZNEW = fc(3,:);
        end


    end
end