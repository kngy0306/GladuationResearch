plotRed = zeros(360, 1);
plotGreen = zeros(360, 1);
plotBlue = zeros(360, 1);

%　桜井さんシミュレーション結果
img_ice5 = imread('./result_Rotation_ice5/iceThick--0.6mm__degree--53.png', 'png'); 
%　氷-空気での反射付け加え結果
img_calc = imread('./result_calculate/iceThick--0.6mm__degree--53.png', 'png');

for i = 1 : 360;
  plotRed(i, 1)    = img_calc(1,  i, 1) - img_ice5(1,  i, 1);
  plotGreen(i, 1) = img_calc(1,  i, 2) - img_ice5(1,  i, 2);
  plotBlue(i, 1)    = img_calc(1,  i, 3) - img_ice5(1,  i, 3);
endfor

clf;
h = bar(plotRed);
set(h, "facecolor", "r");
title('Red');