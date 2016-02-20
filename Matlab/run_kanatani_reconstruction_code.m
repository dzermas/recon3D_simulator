%% Execute the Kanatani C code located in the 3D_reconstruction folder

%% Write the points into .dat files
fileIDw1 = fopen('../kanatani_codes/projective_reconstruction/truck/truck1.dat','w');
fprintf(fileIDw1,'%5i\n',length(v_k));
fprintf(fileIDw1,'%8.3f %8.3f\n',v_k(1,:),v_k(2,:));
fclose(fileIDw1);

fileIDw2 = fopen('../kanatani_codes/projective_reconstruction/truck/truck2.dat','w');
fprintf(fileIDw2,'%5i\n',length(u_k));
fprintf(fileIDw2,'%8.3f %8.3f\n',u_k(1,:),u_k(2,:));
fclose(fileIDw2);

%% Call the C code for the reconstruction!
system('./truck_reconstruction.sh');

%% Read the Projection matrices
fidP1 = fopen('../kanatani_codes/projective_reconstruction/truck/result/proj-000.dat','r');
P1_cell = textscan(fidP1, '%f%f%f%f');
fclose(fidP1);

P1 = [P1_cell{1} P1_cell{2} P1_cell{3} P1_cell{4}];

fidP2 = fopen('../kanatani_codes/projective_reconstruction/truck/result/proj-001.dat','r');
P2_cell = textscan(fidP2, '%f%f%f%f');
fclose(fidP2);

P2 = [P2_cell{1} P2_cell{2} P2_cell{3} P2_cell{4}];

%% Read the 3D points
fid = fopen('../kanatani_codes/projective_reconstruction/truck/result/3D.dat','r');
datacell = textscan(fid, '%f%f%f%f', 'HeaderLines', 1);
fclose(fid);

A = [datacell{1} datacell{2} datacell{3} datacell{4}]';

Y1 = P1*A;
Y2 = P2*A;




for i=1:length(A)
    Z(i) = 1/A(4,i);
end

for i=1:length(A)
    B(1,i) = A(1,i)/A(4,i);
    B(2,i) = A(2,i)/A(4,i);
    B(3,i) = A(3,i)/A(4,i);
end

figure, plot3(Y2(1,:)', Y2(2,:)', Z, 'r*')
 
figure, hold on, grid on
plot3(B(1,:),B(2,:),B(3,:),'g.')
plot3(Y1(1,:)/Y1(3,:),Y1(2,:)/Y1(3,:),Y1(3,:)/Y1(3,:),'r.')
plot3(Y2(1,:),Y2(2,:),Y2(3,:),'b.')
plot3(A(1,:),A(2,:),A(3,:),'g.')

figure, hold on
for j=1:length(Y1)
    psi = ceil(abs(v_k(1,j)));
    chi = ceil(abs(u_k(2,j)));
    plot3(Y1(1,j),Y1(2,j),Z(j),'Marker','*',...
       'Color',[image1(chi,psi,3)/255, image1(chi,psi,2)/255, image1(chi,psi,1)/255])
end
    
    
    