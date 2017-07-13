function [E, residual] = computeEssentialHL(C1, C2)
% Higgins-Longuet or 8 point algorithm for estimating the Fundamental matrix
% which is the Essential Matrix here, since we have normalized the 2D
% keypoints w.r.t. the camera intrinsics.
% Input: 
% C1, C2: 2D corresponding keypoints in the two images, if the points are
% normalized we get back the Essential matrix
% Output
% E: Fundamental/Essential matrix 3x3
% residual: the average residual of the epipolar constraint for all points
% used to compute the fundamental/essential matrix

A = zeros(size(C1,2), 9);
for i=1:size(C1,2)
    u = C1(:,i);
    v = C2(:,i);
    A(i,:) = [u(1)*v(1) u(1)*v(2) u(1) u(2)*v(1) u(2)*v(2) u(2) v(1) v(2) 1];
end

x = solveHomogeneousEq(A);
E = [x(1:3)'; x(4:6)'; x(7:9)'];
% SVD cleanup
[U, D, V] = svd(E);
D(1,1) = 1;
D(2,2) = 1;
D(3,3) = 0;
E = U*D*V';

ss = 0;
for i=1:length(C1)
    ss = ss + C2(:,i)'*E*C1(:,i);
    % assert (u_c(:,i)'*G*v_c(:,i) < 10e-03, 'Error Fundamental');
end
residual = ss / length(C1);