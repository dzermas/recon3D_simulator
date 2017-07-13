function E = ComputeEssentialHL(C1, C2)
% Higgins-Longuet or 8 point algorithm for computing the Fundamental matrix
% which is the Essential Matrix here, since we have normalized the 2D
% keypoints w.r.t. the camera intrinsics.
% Input: 2D corresponding keypoints in the two images C1 and C2
% Output: Fundamental matrix 3x3

A = zeros(size(C1,2), 9);
for i=1:size(C1,2)
    u = C1(:,i);
    v = C2(:,i);
    A(i,:) = [u(1)*v(1) u(2)*v(1) v(1) u(1)*v(2) u(2)*v(2) v(2) u(1) u(2) 1];
end

x = SolveHomogeneousEq(A);
E = [x(1:3)'; x(4:6)'; x(7:9)'];
% SVD cleanup
[U, D, V] = svd(E);
D(1,1) = 1;
D(2,2) = 1;
D(3,3) = 0;
E = U*D*V';
