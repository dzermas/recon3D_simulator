function [E_RANSAC, keys1, keys2] = RansacOutlierRejection(ka, kb)

% each row index in the two keypoint 

% RANSAC threshold
threshold = 10^(-2);
% RANSAC success rate
p = 0.95;
% Probability of choosing an inlier
w = 0.8;
% Number of samples to build the model
n = 8;
% Number of itearations for RANSAC
nIter = log(1-p) / log(1 - w^n);
% Max Inliers

max_nInliers = 0;
inlierIdx = [];
for i=1:nIter
    % Select 8 random points
    randIdx = randperm(size(ka,2), n);
    % Compute Essential Matrix model
    E = ComputeEssentialHL(ka(:,randIdx), kb(:,randIdx));

    % Test Essential model
    nInliers = 0;
    for j = 1 : size(ka,2)
        e = abs([ka(:,j); 1]' * E * [kb(:,j); 1]);
        if (e < threshold)
            nInliers = nInliers + 1;
            inlierIdx = [inlierIdx; j];
        end
    end
    if (nInliers > max_nInliers) 
        max_nInliers = nInliers; 
        E_RANSAC = E;
    end
    
end

keys1 = ka(:,inlierIdx);
keys2 = kb(:,inlierIdx);