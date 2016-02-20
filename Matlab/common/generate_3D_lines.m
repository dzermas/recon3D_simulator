function lines_3D = generate_3D_lines(N_lines, depth_range)
%% Generates the end points of the 3D line segments:
% (x1,y1,z1) -> (x2,y2,z2)

for i=1:N_lines
    lines_3D.x1(i) = depth_range * (2*rand - 1);
    lines_3D.x2(i) = depth_range * (2*rand - 1);
    lines_3D.y1(i) = depth_range * (2*rand - 1);
    lines_3D.y2(i) = depth_range * (2*rand - 1);
    lines_3D.z1(i) = depth_range * 2*rand;
    lines_3D.z2(i) = depth_range * 2*rand;
end