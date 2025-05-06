function [dists, nearestInds] = searchOctTreeRadius(point, radius, points, binParents, binCorners, pointBins)

% Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

% Last modified at 15.36 on 29/04/2025


% Function finds nearest neighbours within a search radius

[pointBin] = findPointBin(point, binParents, binCorners);
searchX = point(1); searchY = point(2); searchZ = point(3);
dists = []; nearestInds = [];
[dists, nearestInds] = newNearCheck(searchX, searchY, searchZ, points, pointBins, dists, nearestInds, pointBin, radius);

checkBins = []; % Bins to be checked
parent = binParents(pointBin); % Parent of start bin
siblings = find(binParents == parent); % Siblings of start bin
siblings(siblings==pointBin) = []; % Removing start bin from siblings list

checkBins(end+1:end+length(siblings)) = siblings;


counter = 1;

% checkBins runs ahead of the actual bin being checked. Once the bin being
% checked has caught up to the checkBins then there are no more bins to
% check and the algorithm automatically stops

while counter <= length(checkBins)

    checkBin = checkBins(counter); % Bin being examined
    dist = point2binDist(searchX, searchY, searchZ, checkBin, binCorners); % Check distance to considered bin

    % If the distance is less than the maximum distance see if bin has children and add to checking list
    if dist < radius
        children = find(binParents == checkBin);
        if any(children)
            checkBins(end+1:end+length(children)) = children;
        else
        % RECOMPUTE NEAREST NEIGHBOUR HERE
        [dists, nearestInds] = newNearCheck(searchX, searchY, searchZ, points, pointBins, dists, nearestInds, checkBin, radius);
        end
    end


    % If the bin does not have a parent of 1, then find its siblings and add to list of bins to check
    if parent ~=1
        oldParent = parent;
        parent = binParents(parent); % Finding parent of current parent
        siblings = find(binParents == parent); % Siblings of new parent bin
        siblings(siblings==oldParent) = []; % Removing old parent bin from siblings list
        checkBins(end+1:end+length(siblings)) = siblings; % Adding the siblings to the checking list
    end

    % Stepping through array of bins to check
    counter = counter + 1;
end
end


%% Functions

% Function to find which bin a point is located in
function [start] = findPointBin(point, binParents, binCorners)
    % Unpacking coordinates 
    searchX = point(1);
    searchY = point(2);
    searchZ = point(3);
    start = 0; % Algorithm starts in bin 0
    while true
        check = find(binParents == start); % Indexes of child bins
        if ~any(check) % If there are no child bins
            break
        end
        for j = 1:length(check) % Looping through all child bins and checking if need to look further
            if searchX >= binCorners(check(j),1) && searchX <=  binCorners(check(j),4) && searchY >=  binCorners(check(j),2) && searchY <=  binCorners(check(j),5) && searchZ >=  binCorners(check(j),3) && searchZ <=  binCorners(check(j),6)
                start = check(j); % Setting child bin to examine again
            end
        end
    end
end


function [dists, nearestInds] = newNearCheck(searchX, searchY, searchZ, points, pointBins, dists, nearestInds, checkBin, radius)
    binpointInds = find(pointBins==checkBin); % All points indexes in the bin
    binpoints = points(binpointInds,:); % All points in the bin
    newdists = sqrt((searchX - binpoints(:,1)).^2 + (searchY - binpoints(:,2)).^2 + (searchZ - binpoints(:,3)).^2); % Dist from target to all points in bin
    sel = newdists > radius;
    newdists(sel) = []; binpointInds(sel) = [];
    dists(end+1:end+length(newdists),1) = newdists;
    nearestInds(end+1:end+length(binpointInds),1) = binpointInds;
end


function [dist] = point2binDist(targX, targY, targZ, examBin, binCorners)

    minX = binCorners(examBin,1);
    minY = binCorners(examBin,2);
    minZ = binCorners(examBin,3);
    maxX = binCorners(examBin,4);
    maxY = binCorners(examBin,5);
    maxZ = binCorners(examBin,6);

     % Computing distance to potential bins to check
        if targX < minX && targY > minY && targY < maxY && targZ >= maxZ % TML
            dist = sqrt((targX-minX)^2 + (targZ-maxZ)^2);
        elseif targX <= minX && targY >= maxY && targZ >= maxZ % TTL
            dist = sqrt((targX-minX)^2 + (targY-maxY)^2 + (targZ-maxZ)^2);
        elseif targX > minX && targX < maxX && targY > maxY && targZ >= maxZ % TTM
            dist = sqrt((targY-maxY)^2 + (targZ-maxZ)^2);
        elseif targX >= maxX && targY >= maxY && targZ >= maxZ % TTR
            dist = sqrt((targX-maxX)^2 + (targY-maxY)^2 + (targZ-maxZ)^2);
        elseif targX > maxX && targY > minY && targY < maxY && targZ >= maxZ % TMR
            dist = sqrt((targX-maxX)^2 + (targZ-maxZ)^2);
        elseif targX >= maxX && targY <= minY && targZ >= maxZ % TBR
            dist = sqrt((targX-maxX)^2 + (targY-minY)^2 + (targZ-maxZ)^2);
        elseif targX > minX && targX < maxX && targY < minY && targZ >= maxZ % TBM
            dist = sqrt((targY-minY)^2 + (targZ-maxZ)^2);
        elseif targX <= minX && targY <= minY && targZ >= maxZ % TBL
            dist = sqrt((targX-minX)^2 + (targY-minY)^2 + (targZ-maxZ)^2);
        elseif targX >= minX && targX <= maxX && targY >= minY && targY <= maxY && targZ >= maxZ % TMM
            dist = maxZ - targZ;

            

        elseif targX <= minX && targY > minY && targY < maxY && targZ > minZ  && targZ < maxZ % MML
            dist = minX - targX;
        elseif targX <= minX && targY >= maxY && targZ > minZ  && targZ < maxZ % MTL
            dist = sqrt((targX-minX)^2 + (targY-maxY)^2);
        elseif targX > minX && targX < maxX && targY >= maxY && targZ > minZ  && targZ < maxZ % MTM
            dist = maxY - targY;
        elseif targX >= maxX && targY >= maxY && targZ > minZ  && targZ < maxZ % MTR
            dist = sqrt((targX-maxX)^2 + (targY-maxY)^2);
        elseif targX >= maxX && targY > minY && targY < maxY && targZ > minZ  && targZ < maxZ % MMR
            dist = maxX - targX;
        elseif targX >= maxX && targY <= minY && targZ > minZ  && targZ < maxZ % MBR
            dist = sqrt((targX-maxX)^2 + (targY-minY)^2);
        elseif targX > minX && targX < maxX && targY <= minY && targZ > minZ  && targZ < maxZ % MBM
            dist = minY - targY;
        elseif targX <= minX && targY <= minY && targZ > minZ  && targZ < maxZ % MBL
            dist = sqrt((targX-minX)^2 + (targY-minY)^2);



        elseif targX < minX && targY > minY && targY < maxY && targZ <= minZ % BML
            dist = sqrt((targX-minX)^2 + (targZ-minZ)^2);
        elseif targX <= minX && targY >= maxY && targZ <= minZ % BTL
            dist = sqrt((targX-minX)^2 + (targY-maxY)^2 + (targZ-minZ)^2);
        elseif targX > minX && targX < maxX && targY > maxY && targZ <= minZ % BTM
            dist = sqrt((targY-maxY)^2 + (targZ-minZ)^2);
        elseif targX >= maxX && targY >= maxY && targZ <= minZ % BTR
            dist = sqrt((targX-maxX)^2 + (targY-maxY)^2 + (targZ-minZ)^2);
        elseif targX > maxX && targY > minY && targY < maxY && targZ <= minZ % BMR
            dist = sqrt((targX-maxX)^2 + (targZ-minZ)^2);
        elseif targX >= maxX && targY <= minY && targZ <= minZ % BBR
            dist = sqrt((targX-maxX)^2 + (targY-minY)^2 + (targZ-minZ)^2);
        elseif targX > minX && targX < maxX && targY < minY && targZ <= minZ % BBM
            dist = sqrt((targY-minY)^2 + (targZ-minZ)^2);
        elseif targX <= minX && targY <= minY && targZ <= minZ % BBL
            dist = sqrt((targX-minX)^2 + (targY-minY)^2 + (targZ-minZ)^2);
        elseif targX >= minX && targX <= maxX && targY >= minY && targY <= maxY && targZ <= minZ % BMM
            dist = minZ - targZ;
        end

dist = abs(dist);
end
