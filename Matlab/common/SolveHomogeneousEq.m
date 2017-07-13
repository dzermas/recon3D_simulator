function x = SolveHomogeneousEq(A)
% The solution of a homogeneous equation Ax = 0, is the nullspace of A, or
% the eigenvector of the smallest eigenvalue.
[~,~,V] = svd(A);
x = V(:,end);