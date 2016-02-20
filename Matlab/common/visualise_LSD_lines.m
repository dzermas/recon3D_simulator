%% Visualize line segments in images
function visualise_LSD_lines(lines, image)

% Plot the lines!
figure, imshow(image), hold on;
for k = 1:size(lines,1)
    dist = sqrt((lines(k,1)-lines(k,3))^2 + (lines(k,2)-lines(k,4))^2);
    if dist > 0
       xy = [lines(k,1) lines(k,2); lines(k,3) lines(k,4)];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
end
