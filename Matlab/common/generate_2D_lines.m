function lines_2D = generate_2D_lines(P, lines_3D)

for i=1:size(lines_3D.x1,2)
     tmp = P * [lines_3D.x1(i) lines_3D.y1(i) lines_3D.z1(i) 1]';
     lines_2D.x1(i) = tmp(1)/tmp(3);
     lines_2D.y1(i) = tmp(2)/tmp(3);
     lines_2D.z1(i) = tmp(3)/tmp(3);
     
     tmp = P * [lines_3D.x2(i) lines_3D.y2(i) lines_3D.z2(i) 1]';
     lines_2D.x2(i) = tmp(1)/tmp(3);
     lines_2D.y2(i) = tmp(2)/tmp(3);
     lines_2D.z2(i) = tmp(3)/tmp(3);
end
