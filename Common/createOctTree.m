function [binParents, binCorners, pointBins] = createOctTree(points, npointsbin, initialBounds)

    % Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
    % of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
    % Energy.
    
    % Copyright under a BSD 3-Clause License, see
    % https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git
    
    % Last modified at 15.38 on 29/04/2025


    % Function to be called by user to create the octree
    
    pointBins = ones(size(points,1),1); % Stores which bin each point belongs to
    binParents = [0]; % Storing where each point came from
    binCorners = initialBounds;
    binNo = 1;
    [binParents, binCorners, pointBins] = OctTreeRecursive(binNo, binParents, binCorners, pointBins, points, npointsbin);
end

% Funtcion called by createOctTree to recursively create the tree
function [binParents, binCorners, pointBins] = OctTreeRecursive(binNo, binParents, binCorners, pointBins, points, npointsbin)

binXl = (binCorners(binNo,4) - binCorners(binNo,1))/2;
binYl = (binCorners(binNo,5) - binCorners(binNo,2))/2;
binZl = (binCorners(binNo,6) - binCorners(binNo,3))/2;
newBinNo = size(binCorners,1); % Updating the newBinNo to the last bin created

for j = 1:8
    if j==1
        binCorner = [binCorners(binNo,1) binCorners(binNo,2) binCorners(binNo,3)+binZl binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,6)];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==2
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2) binCorners(binNo,3)+binZl binCorners(binNo,4) binCorners(binNo,2)+binYl binCorners(binNo,6)];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==3
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2) binCorners(binNo,3) binCorners(binNo,4) binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==4
        binCorner = [binCorners(binNo,1) binCorners(binNo,2) binCorners(binNo,3) binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==5
        binCorner = [binCorners(binNo,1) binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl binCorners(binNo,1)+binXl binCorners(binNo,5) binCorners(binNo,6)];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==6
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl binCorners(binNo,4) binCorners(binNo,5) binCorners(binNo,6)];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==7
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3) binCorners(binNo,4) binCorners(binNo,5) binCorners(binNo,3)+binZl];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==8
        binCorner = [binCorners(binNo,1) binCorners(binNo,2)+binYl binCorners(binNo,3) binCorners(binNo,1)+binXl binCorners(binNo,5) binCorners(binNo,3)+binZl];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    end


        newBinNo = newBinNo + 1; % Increasing the bin number if there is anything in this bin
        binParents(newBinNo,1) = binNo; % Storing where the new bins came from
        binCorners(newBinNo,1:6) = binCorner; % Storing the corners of each bin


        pointBins(sel) = newBinNo; % Updating which subdivision a point lies in

        if length(find(sel)) > npointsbin
            [binParents, binCorners, pointBins] = OctTreeRecursive(newBinNo,binParents, binCorners, pointBins, points, npointsbin);
            newBinNo = size(binCorners,1); % Updating the newBinNo to the last bin created
        end

end
end