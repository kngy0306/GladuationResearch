%豌ｷ縺ｮ繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ
clear
tic; %螳溯｡梧凾髢捺ｸｬ螳�

%--------------------------------繝�く繧ｹ繝医ヵ繧｡繧､繝ｫ縺ｮ隱ｭ縺ｿ霎ｼ縺ｿ----------------------------
%780nm縺九ｉ380nm縺ｾ縺ｧ縺ｮ謨ｰ蛟､縺瑚ｨ伜�縺輔ｌ縺ｦ縺�ｋ繝�く繧ｹ繝医ヵ繧｡繧､繝ｫ
x1=load('x-lambda.txt');
y1=load('y-lambda.txt');
z1=load('z-lambda.txt');
S=load('D65.txt');

%-------------------------------逕ｻ蜒上ヵ繧｡繧､繝ｫ縺ｮ隱ｭ縺ｿ霎ｼ縺ｿ--------------------------------
img=imread('Chromaticity_Diagram_format.png','png'); %xy濶ｲ蠎ｦ蝗ｳ縺ｮ逕ｻ蜒剰ｪｭ縺ｿ霎ｼ縺ｿ
img=double(img)/255;%(0��255縺ｮ謨ｰ蛟､縺ｧ縺ｧ隱ｭ縺ｿ霎ｼ繧薙□縺ｮ繧偵�0.0��1.0縺ｫ螟画鋤)

%----------xy濶ｲ蠎ｦ蝗ｳ縺ｮ菴咏區縺ｮ髟ｷ縺�----------
%菴咏區�嘖pace
sx=0;
sy=0;
for k=1:70% x霆ｸ
    if img(450,k)<0.5 sx=k; end
end
for k=599:-1:530% y霆ｸ
    if img(k,200)<0.5 sy=599-k; end
end

%----------xy濶ｲ蠎ｦ蝗ｳ�醍岼逶帙ｊ縺ｮ髟ｷ縺�----------
%逶ｮ逶帙ｊ�嘖cale
x1s=0; % 蛻晄悄蛟､
for k=55:90
    if img(450,k)<0.9 x1s=k-sx; end
end

xy=zeros(180,2);

%----------蜃ｺ蜉帷判蜒上�----------
height= 50;
width = 360;
result1=zeros(height,width,3); % 邨先棡繧定｡ｨ縺咏判蜒上ｒ蛻晄悄蛹�  (邵ｦ,讓ｪ,RGB)縺ｮ蝗ｳ縺ｮ繧ｵ繧､繧ｺ
result2=zeros(height,width,3); 

%---------------蛻晄悄險ｭ螳�------------------------------------------------
% 蝗櫁ｻ｢  �� rotation[rp1��,rb1�枉
% 蛛丞�譚ｿ�� polarizer[P1�枉
% 隍�ｱ域釜�� birefringence[Ba1�枉

%0.005

%fname1 = '豌ｷ縺ｮ繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ邨先棡.png';
%fname2 = '驕�嶌霆ｸ蛯ｾ譁應ｸｭ蠢�ｻｸ.png';
%xyfname1 = '譎ｮ騾壹�蝗櫁ｻ｢濶ｲ蠎ｦ蝗ｳ845.png';
%xyfname2 = '蝗櫁ｻ｢濶ｲ蠎ｦ蝗ｳ845(10ﾂｰ縺斐→).png';

%-------------蛛丞�譚ｿ縺ｮ險ｭ螳�---------------
rP = [0 pi/2 pi 0 0 0 0 0 0 0];% 蛛丞�譚ｿ縺ｮ隗貞ｺｦ(10譫壼�
p=[cos(rP(1))^2 sin(rP(1))*cos(rP(1));sin(rP(1))*cos(rP(1)) sin(rP(1))^2];% 蛛丞�譚ｿ縺ｮ陦悟�險育ｮ�
p2=[cos(rP(2))^2 sin(rP(2))*cos(rP(2));sin(rP(2))*cos(rP(2)) sin(rP(2))^2];% 蛛丞�譚ｿ縺ｮ陦悟�險育ｮ�

%-------------菴咲嶌蟾ｮ繝輔ぅ繝ｫ繝�縺ｮ險ｭ螳�-----------------
rB=[pi/4 pi*3/4 0 0 0 0 0 0 0 0];% 隍�ｱ域釜諤ｧ譚先侭縺ｮ隗貞ｺｦ(10譫�)

Ba1=[cos(rB(1)) -sin(rB(1));sin(rB(1)) cos(rB(1))];% 隍�ｱ域釜諤ｧ譚先侭縺ｮ蝗櫁ｻ｢險育ｮ�
Ba2=[cos(rB(2)) -sin(rB(2));sin(rB(2)) cos(rB(2))];% 隍�ｱ域釜諤ｧ譚先侭縺ｮ蝗櫁ｻ｢險育ｮ�
Ba3=[cos(rB(3)) -sin(rB(3));sin(rB(3)) cos(rB(3))];% 隍�ｱ域釜諤ｧ譚先侭縺ｮ蝗櫁ｻ｢險育ｮ�


%----------蟷ｲ貂芽牡繝√Ε繝ｼ繝医�濶ｲ蠎ｦ蝗ｳ菴懈�----------
%(1,1��3)縺ｫ豁｣髱｢縺ｮ濶ｲ縲�(2,1��3)縺ｫ譁懊ａ縺ｮ濶ｲ繧剃ｿ晏ｭ倥☆繧�
%(1��2,1)縺ｫR,(1��2,2)縺ｫG,(1��2,3)縺ｫB縺ｮ蛟､繧剃ｿ晏ｭ倥☆繧�

%**********險育ｮ�**********

% 繝ｪ繧ｿ繝��繧ｷ繝ｧ繝ｳ�嗷etardation
% 豕｢髟ｷ�嗹avelength[w]
% 蛻�ｯ搾ｼ單enominator[kd]
% 蛻�ｭ撰ｼ嗜umerator[Xn�杙n]
% 蠎ｧ讓呻ｼ喞oordinate

%----------隕ｳ貂ｬ譚｡莉ｶ------------------%
angle_a = 53;         %蜈･蟆�ｧ�(繝悶Μ繝･繝ｼ繧ｹ繧ｿ繝ｼ隗�:53.1ﾂｰ)
ice_thick = 600000; %豌ｷ縺ｮ蜴壹＆[nm]
%蜈･蟆�婿蜷�(螻域釜邇㌢->螻域釜邇㍍)
%豌ｷ縺ｮ螻域釜邇�(1.309,1.313)
n_air = 1;              %遨ｺ豌励�螻域釜邇�(1)
n_x = 1.309;            %螻域釜邇�1
n_y = 1.313;            %螻域釜邇�2
n_z = 1.309;            %z霆ｸ 蜴壹＆譁ｹ蜷大ｱ域釜邇�(c霆ｸ譁ｹ蜷�)(1.313)
n_water = 1.333;        %豌ｴ縺ｮ螻域釜邇�(1.333)

%----------險育ｮ�----------------------%
rad_a = deg2rad(angle_a);                   %蜈･蟆�ｧ�(繝ｩ繧ｸ繧｢繝ｳ)
rad_b = asin( sin(rad_a) * (n_air/n_y) );     %螻域釜隗�(繝ｩ繧ｸ繧｢繝ｳ)
angle_b = rad2deg(rad_b);   %蜃ｺ蟆�ｧ�(螻域釜隗�)
n_x = n_x * n_z / (sqrt(n_x.^2 * sin(rad_b).^2 + n_z.^2 * cos(rad_b).^2)); %蜈･蟆�ｧ貞ｱ域釜邇�
n_y = n_y * n_z / (sqrt(n_y.^2 * sin(rad_b).^2 + n_z.^2 * cos(rad_b).^2)); %蜈･蟆�ｧ貞ｱ域釜邇�

rs = Rs(angle_a,n_air,n_water);
rp = Rp(angle_a,n_air,n_water);
percentage = rs / (rs+rp);   %s蛛丞�縺ｮ蜑ｲ蜷�
%----------繝ｪ繧ｿ驥� 竊� xyz 竊� xy,rgb ----------

% 縺薙＞縺､繧峨�險育ｮ励〒荳譎ら噪縺ｫ菴ｿ縺�ユ繝ｳ繝昴Λ繝ｪ繝ｼ螟画焚
% 險育ｮ礼ｵ先棡縺ｯ陦悟�縺ｫ菫晏ｭ倥☆繧�

%r = [140 565 705 845];% 髢句ｧ九Μ繧ｿ繝��繧ｷ繝ｧ繝ｳ驥充nm]..............菴咲嶌蟾ｮ繝輔ぅ繝ｫ繝�縺ｮ遞ｮ鬘�(1��3縺ｧ驕ｸ縺ｹ繧�)
ice = abs((n_y - n_x) * (ice_thick / cos(rad_b)));
%ice = 0.00314 * (ice_thick/0.885528); %譁懊ａ37.5ﾂｰ縺九ｉ蜈･蟆�凾繝ｪ繧ｿ驥充nm] 騾ｲ螻郁ｻｸ
%ice2 = 0.004 * (ice_thick/0.885528); %驕�ｱ郁ｻｸ
%ice3 = 0.00086 * (ice_thick/0.885528); %蜴壹＆譁ｹ蜷�

% X�杙 蛛丞�迥ｶ諷� 蜃ｺ蟆��蠑ｷ蠎ｦ繧貞�譛溷喧
X  = 0;
Y  = 0;
Z  = 0;
X2 = 0;
Y2 = 0;
Z2 = 0;
kd = 0;
I  = zeros(81,1);
I2 = zeros(81,1);

for ri = 1:360;
    
    X  = 0;
    Y  = 0;
    Z  = 0;
    X2 = 0;
    Y2 = 0;
    Z2 = 0;
    
    kd = 0;
    
    for t = 1:81% MATLAB縺ｮ豺ｻ縺亥ｭ励�1縺九ｉ縺ｪ縺ｮ縺ｧ縲√％縺ｮ繧医≧縺ｫ縺励※縺�ｋ
        w = t*5+375;%豕｢髟ｷ380��780nm縺ｾ縺ｧ(蜿ｯ隕門�繧�5[nm]縺壹▽隱ｿ縺ｹ繧�)
        
        %F1=[exp(-1i*pi*r(1)/w) 0;0 exp(1i*pi*r(1)/w)];% ﾎｻ譚ｿ縺ｮ豕｢髟ｷ豈弱�陦悟�險育ｮ�
        %F2=[exp(-1i*pi*r(2)/w) 0;0 exp(1i*pi*r(2)/w)];% ﾎｻ/4譚ｿ縺ｮ豕｢髟ｷ豈弱�陦悟�險育ｮ�
        F_ice=[exp(-1i*pi*ice/w) 0;0 exp(1i*pi*ice/w)]; %豌ｷ縺ｮ豕｢髟ｷ縺斐→縺ｮ陦悟�險育ｮ�
        %F_ice2=[exp(-1i*pi*ice2/w) 0;0 exp(1i*pi*ice2/w)];
        %F_ice3=[exp(-1i*pi*ice3/w) 0;0 exp(1i*pi*ice3/w)];
        P=[cos(pi*(ri-1)/180)^2 sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180);sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180) sin(pi*(ri-1)/180)^2];% 蛛丞�譚ｿ縺ｮ陦悟�險育ｮ�
        
        %a�柬縺ｾ縺ｧ菴ｿ縺��繧貞酔縺倥ｈ縺�↓譖ｸ縺�
        %E=P1 *Ba1*F*Ba2 *Bb1*F4*Bb2 *Bc1*F*Bc2 *Bd1*F4*Bd2 *P2; %蜃ｺ蟆��蠑ｷ蠎ｦ縺ｮ險育ｮ�,繧ｸ繝ｧ繝ｼ繝ｳ繧ｺ繝槭ヨ繝ｪ繝�け繧ｹ
        Es = P * (Ba1*F_ice*Ba1') * [1;0];
        Ep = P * (Ba1*F_ice*Ba1') * [0;1];
        %E  = p  *(Ba1*F1*Ba1')*P;
        %E2 = p2 *(Ba1*F1*Ba1')*P;
        %E3 = p  *(Ba1*F2*Ba1')*P;
        %E4 = p2 *(Ba1*F2*Ba1')*P;
        I(t,1) =sum(percentage * abs(Es(:).^2) + (1-percentage) * abs(Ep(:).^2));
        %I2(t,1)=sum(percentage * abs(E3(:).^2)+(1-percentage) * abs(E4(:).^2));
        
        X=X+(S(t,1)*I(t,1)*x1(t,1));
        Y=Y+(S(t,1)*I(t,1)*y1(t,1)); % Y numerator 蛻�ｭ舌�遨榊�縺ｮ險育ｮ�
        Z=Z+(S(t,1)*I(t,1)*z1(t,1)); % Z numerator 蛻�ｭ舌�遨榊�縺ｮ險育ｮ�
        
        %X2=X2+(S(t,1)*I2(t,1)*x1(t,1));
        %Y2=Y2+(S(t,1)*I2(t,1)*y1(t,1)); % Y numerator 蛻�ｭ舌�遨榊�縺ｮ險育ｮ�
        %Z2=Z2+(S(t,1)*I2(t,1)*z1(t,1)); % Z numerator 蛻�ｭ舌�遨榊�縺ｮ險育ｮ�
        
        kd=kd+(S(t,1)*y1(t,1));% k denominator 蛻�ｯ阪�遨榊�縺ｮ險育ｮ� 100縺ｯ蠕後〒豸医∴繧九�縺ｧ辟｡隕�

    end

    X = X/kd;
    Y = Y/kd;
    Z = Z/kd;
    %X2 = X2/kd;
    %Y2 = Y2/kd;
    %Z2 = Z2/kd;
    
    %kd縺ｧ蜑ｲ縺｣縺ｦ菫晏ｭ�
    xy(ri,1)=X/(X+Y+Z);% x蠎ｧ讓�
    xy(ri,2)=Y/(X+Y+Z);% y蠎ｧ讓�
    %xy2(ri,1)=X2/(X2+Y2+Z2);% x蠎ｧ讓�
    %xy2(ri,2)=Y2/(X2+Y2+Z2);% y蠎ｧ讓�
    
    %----------XYZ蠎ｧ讓吶ｒRGB螟画鋤----------
    % CIE XYZ邉ｻ縺九ｉCIE RGB邉ｻ縺ｸ縺ｮ螟画鋤
    R =    2.3655*X +  (-0.8971)*Y +  (-0.4683)*Z;
    G = (-0.5151)*X +     1.4264*Y +     0.0887*Z;
    B =    0.0052*X +  (-0.0144)*Y +     1.0089*Z;
    %R2 =    2.3655*X2 +  (-0.8971)*Y2 +  (-0.4683)*Z2;
    %G2 = (-0.5151)*X2 +     1.4264*Y2 +     0.0887*Z2;
    %B2 =    0.0052*X2 +  (-0.0144)*Y2 +     1.0089*Z2;
    
    R = max(R,0.0);
    G = max(G,0.0);
    B = max(B,0.0);
    %R2 = max(R2,0.0);
    %G2 = max(G2,0.0);
    %B2 = max(B2,0.0);

    %----------蟷ｲ貂芽牡----------
    result1(:,ri,1)=R; % R縺ｮ蛟､
    result1(:,ri,2)=G; % G縺ｮ蛟､
    result1(:,ri,3)=B;% B縺ｮ蛟､
    %result2(:,ri,1)=R2; % R縺ｮ蛟､
    %result2(:,ri,2)=G2; % G縺ｮ蛟､
    %result2(:,ri,3)=B2; % B縺ｮ蛟､
    
    %----------濶ｲ蠎ｦ蝗ｳ-----------------
    %豎ゅａ縺溯牡蠎ｦ蠎ｧ讓吶ｒ縲』y濶ｲ蠎ｦ蝗ｳ縺ｮ逕ｻ蜒上�逶ｮ逶帙ｊ縺ｮ髟ｷ縺輔↓蜷医ｏ縺帙ｋ縺溘ａ縺ｮ險育ｮ励ｒ陦後＞縲∬｡悟�xy縺ｫ豎ゅａ縺溷､繧呈�ｼ邏�
    %xy=xy*(x1s*2)*10;
    %xy(ri,1)=xy(ri,1)+double(sx); % 
    %xy(ri,2)=xy(ri,2)+double(sy); % 
    %xy=round(xy);
    %xy2=xy2*(x1s*2)*10;
    %xy2(ri,1)=xy2(ri,1)+double(sx); % 
    %xy2(ri,2)=xy2(ri,2)+double(sy); % 
    %xy2=round(xy2);
    %xy3=xy3*(x1s*2)*10;
    %xy3(ri,1)=xy3(ri,1)+double(sx); % 
    %xy3(ri,2)=xy3(ri,2)+double(sy); % 
    %xy3=round(xy3);
    %xy4=xy4*(x1s*2)*10;
    %xy4(ri,1)=xy4(ri,1)+double(sx); % 
    %xy4(ri,2)=xy4(ri,2)+double(sy); % 
    %xy4=round(xy4);

    %599-xy(ri,2)
    %xy濶ｲ蠎ｦ蝗ｳ縺ｫ蛟､繧偵�繝ｭ繝�ヨ
    %img(599-xy(ri,2), xy(ri,1), 1:3)=0;%濶ｲ蠎ｦ蠎ｧ讓吩ｸ翫↓蛟､繧偵�繝ｭ繝�ヨ
    %img2(599-xy2(ri,2), xy2(ri,1), 1:3)=0;
end

%----------逕ｻ蜒丈ｿ晏ｭ�----------
%imwrite(result1,['C:/Users/New User/Desktop/[譯應ｺ評遐皮ｩｶ雉�侭/M2/繧ｸ繝･繧ｨ繝ｪ繝ｼ繝舌ヶ繝ｫ/繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ邨先棡(Rotation_ice4.m)/蜴壹＆',num2str(ice_thick/1000000),'mm_蜈･蟆�ｧ�',num2str(angle_a),'ﾂｰ.png']);
%imwrite(result2,fname2);
%imwrite(img,xyfname1);
imwrite(result1,['./',num2str(ice_thick/1000000),'mm_蜈･蟆�ｧ�',num2str(angle_a),'ﾂｰ.png']);

%winopen(['C:/Users/New User/Desktop/[譯應ｺ評遐皮ｩｶ雉�侭/M2/繧ｸ繝･繧ｨ繝ｪ繝ｼ繝舌ヶ繝ｫ/繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ邨先棡(Rotation_ice4.m)/蜴壹＆',num2str(ice_thick/1000000),'mm_蜈･蟆�ｧ�',num2str(angle_a),'ﾂｰ.png']);
winopen(['./',num2str(ice_thick/1000000),'mm_蜈･蟆�ｧ�',num2str(angle_a),'ﾂｰ.png']);
%image(result1);
%winopen(fname1);
%winopen(fname2);
%winopen(xyfname1);

toc; %螳溯｡梧凾髢捺ｸｬ螳�
%**********END**********
