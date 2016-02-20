function [v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani(C1_h, C2_h, E_estimated)
%% Variables
% Returns the corrected point values and the Dx that is needed for other
% corrections (depth etc.)

%% Correct points to accurately validate the epipolar constraint-----------
% Covariances show the variance of the points
% Assume independent noise in each direction on the image (x, y)
k = [0 0 1]';
Vv = eye(3) - k*k';
Vu = eye(3) - k*k';

% Iterate over all the points
% v_c: points in image 1
% u_c: points in image 2
J = []; v_c = []; u_c = [];
for i=1:length(C1_h)
    v = C1_h(:,i); u = C2_h(:,i);
    
    % Compute the error of the points in the two images
    Dv = (u'*E_estimated'*v * Vv * E_estimated * u) / (norm(Vv * E_estimated' * v, 2)^2 + norm(Vu * E_estimated * u,2)^2);
    Du = (u'*E_estimated'*v * Vu * E_estimated' * v) / (norm(Vv * E_estimated' * v, 2)^2 + norm(Vu * E_estimated * u,2)^2);
    
    % Compute the residual to see whether the points correspond or not
    v_c(1:3,i) = v - Dv; u_c(1:3,i) = u - Du;
    Dv_c(1:3,i) = Dv; Du_c(1:3,i) = Du;
    J(i) = (v'*E_estimated*u)^2 / (norm(Vv * E_estimated' * v_c(1:3,i), 2)^2 + norm(Vu * E_estimated * u_c(1:3,i),2)^2);
end


%% Is kept with confidence 95% if J < CritVal (Degrees of Freedom = 1) ----
keepers = J < 0.04;
sum(keepers);

v_k = v_c(:, keepers); u_k = u_c(:, keepers);
Dv_k = Dv_c(:, keepers); Du_k = Du_c(:, keepers);