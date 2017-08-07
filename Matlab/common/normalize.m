function x_n = normalize(x_kk, camera_params)

%Computes the normalized coordinates xn given the pixel coordinates x_kk
%and the intrinsic camera parameters fc, cc and kc.
%
%INPUT: x_kk: Feature locations on the images
%       camera_params: Structure holding the camera parameters:
%       camera_params.fc: Camera focal length
%       camera_params.cc: Principal point coordinates
%       camera_params.kc: Distortion coefficients
%       camera_params.alpha_c: Skew coefficient
%
%OUTPUT: x_n: Normalized feature locations on the image plane (a 2XN matrix)
%
%Important functions called within that program:
%
%comp_distortion_oulu: undistort pixel coordinates.

% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - camera_params.cc(1))/camera_params.fc(1);(x_kk(2,:) - camera_params.cc(2))/camera_params.fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - camera_params.alpha_c * x_distort(2,:);

if norm(camera_params.kc) ~= 0
	% Third: Compensate for lens distortion:
	x_n = comp_distortion_oulu(x_distort,camera_params.kc);
else
   x_n = x_distort;
end
