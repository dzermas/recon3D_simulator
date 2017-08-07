function visualizeCameraPose(camera, R, t, color)

% camera : camera intrinsics structure
% R      : rotation of the camera
% t      : translation of the camera
% color  : 3 component vector (red, green, blue) range [0,1]

fc = camera.fc(1);
w = camera.cc(1)/fc;
h = camera.cc(2)/fc;
f = 1;

c_center = t;

cam_face = [R t] * [-w  w  w -w -w;
                     h  h -h -h  h;
                     f  f  f  f  f;
                     1  1  1  1  1];

cam_base = [c_center(1) cam_face(1,1) c_center(1) cam_face(1,2) c_center(1) cam_face(1,3) c_center(1) cam_face(1,4);
            c_center(2) cam_face(2,1) c_center(2) cam_face(2,2) c_center(2) cam_face(2,3) c_center(2) cam_face(2,4);
            c_center(3) cam_face(3,1) c_center(3) cam_face(3,2) c_center(3) cam_face(3,3) c_center(3) cam_face(3,4)];

plot3([cam_base(1,:) cam_face(1,:)], ...
      [cam_base(2,:) cam_face(2,:)], ...
      [cam_base(3,:) cam_face(3,:)], 'Color', color)