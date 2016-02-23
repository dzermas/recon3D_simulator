function J = jacobianPerspective(p)
% Jacobian is 2x3 since last row is all zeros
x = p(1); y = p(2); z = p(3);

J = [1/z  0  -x/z^2;
      0  1/z -y/z^2];