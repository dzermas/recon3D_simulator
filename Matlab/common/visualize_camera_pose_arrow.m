function visualize_camera_pose_arrow(R, t, color)

% R      : rotation of the camera
% t      : translation of the camera
% color  : 3 component vector (red, green, blue) range [0,1]

c_center = t;

cam_face = [R t] * [0; 0; 1; 1];

plot3([c_center(1) cam_face(1)], ...
      [c_center(2) cam_face(2)], ...
      [c_center(3) cam_face(3)], 'Color', color)