%% Call ASIFT algorithm from the console
% For more info check README.txt in the ASIFT demo folder
% This function is an altered version of the demo_ASIFT from the original
% authors.

function matches = ASIFT(image1, image2)

imgOutVert = 'imgOutVert.png';
imgOutHori = 'imgOutHori.png';
matchings = 'matchings.txt';
keys1 = 'keys1.txt';
keys2 = 'keys2.txt';

% convert the image to png format 
file_img1_png = 'tmpASIFTinput1.png';
file_img2_png = 'tmpASIFTinput2.png';
imwrite(image1, file_img1_png, 'png');
imwrite(image2, file_img2_png, 'png');

% ASIFT command
command_ASIFT = '../demo_ASIFT_src/demo_ASIFT'; 
command_ASIFT = [command_ASIFT ' ' file_img1_png ' ' file_img2_png ' ' ...
  imgOutVert ' ' imgOutHori ' ' matchings ' ' keys1 ' ' keys2];
% 1 for fast normalized points, 0 for slow accurate points
command_ASIFT = [command_ASIFT ' 0'];
	
% get the number of processors 
[s, w] = unix('grep processor /proc/cpuinfo | wc -l');
num_CPUs = str2num(w);
% set the maximum OpenMP threads to the number of processors 
set_threads = sprintf('export OMP_NUM_THREADS=%d;', num_CPUs);
      
command = [set_threads ' ' command_ASIFT];

% Execute command
status = system(command);

if status == 0
    fid = fopen(matchings,'r');
    matches_cell = textscan(fid, '%f%f%f%f');
    fclose(fid);

    matches = [matches_cell{1} matches_cell{2} matches_cell{3} matches_cell{4}];
else
    disp('Error executing ASIFT');
    matches = 0;
end