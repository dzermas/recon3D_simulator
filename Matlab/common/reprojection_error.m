% Computes the reprojection error sum(sum(T*r - p)^2) for 2 cameras
function error = reprojection_error(R1, h1, R2, h2, r, p1, p2)

% R1, R2: Rotation matrices for cameras 1 and 2
% h1, h2: translation vectors for cameras 1 and 2
% r: estimated 3D points
% p1, p2: image plane bearing points as extracted from the feature
% extraction algorithm
%
% error: the squared error

err_sum = 0;
for i=1:size(p1,2)
    cam1 = [R1 h1]*[r(:,i); 1]/norm([R1 h1]*[r(:,i); 1]) - p1(:,i);
    cam2 = [R2' -R2'*h2]*[r(:,i); 1]/norm([R2' -R2'*h2]*[r(:,i); 1]) - p2(:,i);
    err_sum = err_sum + (cam1+cam2);
end

error = norm(err_sum);