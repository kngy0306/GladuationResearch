%氷のシミュレーション
clear
tic;

% 屈折率、ブリュースター角定義
nAir = 1.000;
nIce = 1.309;
nWater = 1.333;
angleBrewster = 53;

%%
% 氷への入射角(ブリュースター角)
% angle_a2i
%%
angle_a2i = 53;

%%
% 氷中の屈折角and入射角
% angle_i2w
%%
angle_i2a = asin(sin(deg2rad(angleBrewster)) * (1.000/1.309));
angle_i2a = rad2deg(angle_i2a); % 37.598


% 空気-氷 の入射光に対する偏光状態
% ※S偏光のままなので着色なし
%rs1 = Rs(angle_a2i, nAir, nIce);
%rp1 = Rp(angle_a2i, nAir, nIce);

% 氷-空気　の入射光に対する偏光状態
rs2 = Rs(angle_i2a, nIce, nAir);
rp2 = Rp(angle_i2a, nIce, nAir);

% 空気-水 入射光に対する偏光状態
rs3 = Rs(angleBrewster, nAir, nWater);
rp3 = Rp(angleBrewster, nAir, nWater);

%----------------テキストファイルの読み込み-----------
x1 = load('x-lambda.txt'); % XYZ表色系等色関数
y1 = load('y-lambda.txt'); 
z1 = load('z-lambda.txt');
S  = load('D65.txt');

%-------------------画像ファイルの読み込み-------------------------
img = imread('Chromaticity_Diagram_format.png', 'png'); 
img = double(img) / 255;

%-----------出力画像----------
height = 50;
width = 360;
result1 = zeros(height, width, 3);
result2 = zeros(height, width, 3); 
%---------------�����ݒ�------------------------------------------------
%---------------初期設定--------------------------
% 複屈折： birefringence[Ba1～]
%-------------�����ܐ��ޗ�(�X)�̐ݒ�-----------------
%-------------複屈折性材料(氷)の設定-----------------
rB = [pi/4 pi*3/4 0 0 0 0 0 0 0 0];
Ba1 = [cos(rB(1)) -sin(rB(1)); sin(rB(1)) cos(rB(1))]; % 複屈折性材料の回転計算

%---------------観測条件----------------
angle_a = 53; % 入射角(ブリュースター角:53.1°)
ice_thick = 800000; % 氷の厚さ[nm]

% 入射方向(屈折率a->屈折率b)
% 氷の屈折率(1.309,1.313)

n_air = 1; % 空気の屈折率(1)
n_x = 1.309;  %屈折率1
n_y = 1.313;  % 屈折率2
n_z = 1.309;  %z軸 厚さ方向屈折率(c軸方向)(1.313)
n_water = 1.333; % 水の屈折率(1.333)

%---------計算-----------
rad_a = deg2rad(angle_a); % 入射角(ラジアン)
rad_b = asin(sin(rad_a) * (n_air/n_y)); % 屈折角(ラジアン) 論文(3.3, 3.5)



angle_b = rad2deg(rad_b); % 出射角(屈折角)


n_y = n_y * n_z / (sqrt(n_y.^2 * sin(rad_b).^2 + n_z.^2 * cos(rad_b).^2)); % 入射角屈折率
n_x = n_x * n_z / (sqrt(n_x.^2 * sin(rad_b).^2 + n_z.^2 * cos(rad_b).^2)); % 入射角屈折率

rs = Rs(angle_a, n_air, n_water); % (3.1)
rp = Rp(angle_a, n_air, n_water); % (3.2)

%percentage1 = rs1 / (rs1 + rp1);
percentage2 = rs2 / (rs2 + rp2); % 氷 - 空気
percentage3 = rs3 / (rs3 + rp3); % 空気 - 水

%percentage = (percentage1 + percentage2 + percentage3) / 3;

%----------リタ量 → xyz → xy,rgb ----------

% こいつらは計算で一時的に使うテンポラリー変数
% 計算結果は行列に保存する

ice = abs((n_y - n_x) * (ice_thick / cos(rad_b))); % 光路差R 論文(3.4)

% X～Z 偏光状態 出射光強度を初期化
X = 0;
X2 = 0;
X3 = 0;
Y = 0;
Y2 = 0;
Y3 = 0;
Z = 0;
Z2 = 0;
Z3 = 0;
kd = 0;
T2 = zeros(81, 1);
T3 = zeros(81, 1);

for ri = 1 : 360;
    
    X2 = 0;
    Y2 = 0;
    Z2 = 0;
    X3 = 0;
    Y3 = 0;
    Z3 = 0;
    kd = 0;
    
    for t = 1 : 81
        w = t * 5 + 375; % 波長380～780nmまで(可視光を5[nm]ずつ調べる)
        
        F_ice = [exp(-1i*pi*ice/w) 0; 0 exp(1i*pi*ice/w)]; %氷の波長ごとの行列計算
        P = [cos(pi*(ri-1)/180)^2 sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180); sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180) sin(pi*(ri-1)/180)^2]; % 偏光板を0°～359°回転
        
        Ep = P * (Ba1*F_ice*Ba1') * [0; 1]; % 偏光状態(3.10)
        Es = P * (Ba1*F_ice*Ba1') * [1; 0]; % 偏光状態(3.11)

        T2(t,1) = sum(percentage2 * abs(Es(:).^2) + (1-percentage2) * abs(Ep(:).^2));
        T3(t,1) = sum(percentage3 * abs(Es(:).^2) + (1-percentage3) * abs(Ep(:).^2));
        
        X2 = X2 + (S(t,1) * T2(t,1) * x1(t,1)); % 分子の積分計算
        Y2 = Y2 + (S(t,1) * T2(t,1) * y1(t,1));
        Z2 = Z2 + (S(t,1) * T2(t,1) * z1(t,1));
        X3 = X3 + (S(t,1) * T3(t,1) * x1(t,1)); % 分子の積分計算
        Y3 = Y3 + (S(t,1) * T3(t,1) * y1(t,1));
        Z3 = Z3 + (S(t,1) * T3(t,1) * z1(t,1));
        
        kd = kd + (S(t, 1) * y1(t, 1)); % 分母の積分の計算 100は後で消えるので無視

    end

    X2 = X2 / kd;
    Y2 = Y2 / kd;
    Z2 = Z2 / kd;
    X3 = X3 / kd;
    Y3 = Y3 / kd;
    Z3 = Z3 / kd;
    
    %X = (X2 + X3) / 2;
    %Y = (Y2 + Y3) / 2;
    %Z = (Z2 + Z3) / 2;
    
    %---------XYZ座標をRGB変換----------
    % CIE XYZ系からCIE RGB系への変換
    R   =     2.3655*X2 +  (-0.8971)*Y2 +  (-0.4683)*Z2;
    G   = (-0.5151)*X2 +      1.4264*Y2 +      0.0887*Z2;
    B   =     0.0052*X2 +  (-0.0144)*Y2 +      1.0089*Z2;
    R1 =     2.3655*X3 +  (-0.8971)*Y3 +  (-0.4683)*Z3;
    G1 = (-0.5151)*X3 +      1.4264*Y3 +      0.0887*Z3;
    B1 =     0.0052*X3 +  (-0.0144)*Y3 +      1.0089*Z3;
    
    R = (R + R1) / 2;
    G = (G + G1) / 2;
    B = (B + B1) / 2;
    %R = (R + R1);
    %G = (G + G1);
    %B = (B + B1);
    
    R = max(R, 0.0);
    G = max(G, 0.0);
    B = max(B, 0.0);

    %---------干渉色----------
    result1(:, ri, 1) = R;
    result1(:, ri, 2) = G;
    result1(:, ri, 3) = B;

end

%画像保存
imwrite(result1, ['./result_calculate/NEW-iceThick--', num2str(ice_thick / 1000000), 'mm__degree--', num2str(angle_a), '.png']);
imshow(result1);
toc;