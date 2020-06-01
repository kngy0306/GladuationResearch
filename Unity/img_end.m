clear
tic;
img_rgb = imread('img.png', 'png'); 
img_input = img_rgb / 255;

height  = 1;
width   = 97200;
img_output = zeros(height,width);

% input       1 x 90 x 3
% output     1 x 97200 

column = 0;
row = 0;
color = 0;
for i = 0 : 97199;
  temporary = rem(i , 270);
  column = idivide(temporary, 3);      %== i / 3
  row  = idivide(i, 270);   %== i / 270
  color = rem(i, 3);      %== i % 3
  
  img_output(1, i+1) = img_input(column+1, row+1, color+1);
end

save('image_end.txt', 'img_output');
toc;