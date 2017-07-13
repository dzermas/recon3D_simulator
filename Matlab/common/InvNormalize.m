function x = invNormalize(x_n, camera_params)

% express points from the image plain coordinates to pixels

% First check if there is lens distortion and add it back.
if norm(camera_params.kc) ~= 0
	% Third: Compensate for lens distortion:
	xd = inv_comp_distortion_oulu(x_n, camera_params.kc);
else
    xd = x_n;
end

% Second apply skew
xd(1,:) = xd(1,:) + camera_params.alpha_c * xd(2,:);

% Third multiply focal length and add principal point
x = [xd(1,:) * camera_params.fc(1) + camera_params.cc(1); xd(2,:) * camera_params.fc(2) + camera_params.cc(2)];