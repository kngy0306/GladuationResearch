clear
tic;

key = 0;
while ( (key < 1) || (key > 10) ) 
  prompt = "Select images 0. ? mm (1~10)\n";
  key = input(prompt);
end


% 1 -> 0.1mmの画像
imgArray = containers.Map({1,2,3,4,5,6,7,8,9,10},
{"ice0_1.png","ice0_2.png","ice0_3.png","ice0_4.png","ice0_5.png","ice0_6.png","ice0_7.png","ice0_8.png","ice0_9.png","ice1_0.png",});

% 1~10mmの画像をimgファイルから取り出す
path = strcat("./img/", imgArray(key));

img_rgb = imread(path, 'png'); 
%img_input = img_rgb / 255; % Octaveで作成した画像の場合
img_input = img_rgb; % MATLABで作成した画像の場合

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

saveName = strcat(num2str(key), '.txt');
outputPath = strcat("./arrayResult/", saveName);
save(outputPath, 'img_output');
toc;