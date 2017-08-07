function C2_R_C1 = generateRandomRotation()

random_angle_1 = rand * 5 * pi / 180;
random_angle_2 = rand * (-5) * pi / 180;
random_angle_3 = rand * 5 * pi / 180;

R_Z = [cos(random_angle_3) -sin(random_angle_3) 0;
    sin(random_angle_3) cos(random_angle_3) 0;
    0 0 1];

R_Y = [cos(random_angle_2) 0 sin(random_angle_2);
        0 1 0;
        -sin(random_angle_2) 0 cos(random_angle_2)];

R_X = [1 0 0;
    0 cos(random_angle_1) -sin(random_angle_1);
    0 sin(random_angle_1) cos(random_angle_1)];

C2_R_C1 = R_Z * R_Y * R_X;