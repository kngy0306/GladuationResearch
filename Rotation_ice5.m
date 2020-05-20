%氷のシミュレーション
clear
tic;

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
ice_thick = 600000; % 氷の厚さ[nm]

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
percentage = rs / (rs + rp); % s偏光の割合
%----------リタ量 → xyz → xy,rgb ----------

% こいつらは計算で一時的に使うテンポラリー変数
% 計算結果は行列に保存する

ice = abs((n_y - n_x) * (ice_thick / cos(rad_b))); % 光路差R 論文(3.4)

% X～Z 偏光状態 出射光強度を初期化
X = 0;
Y = 0;
Z = 0;
kd = 0;
T = zeros(81, 1);

for ri = 1 : 360;
    
    X = 0;
    Y = 0;
    Z = 0;
    kd = 0;
    
    for t = 1 : 81
        w = t * 5 + 375; % 波長380～780nmまで(可視光を5[nm]ずつ調べる)
        
        F_ice = [exp(-1i*pi*ice/w) 0; 0 exp(1i*pi*ice/w)]; %氷の波長ごとの行列計算
        P = [cos(pi*(ri-1)/180)^2 sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180); sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180) sin(pi*(ri-1)/180)^2]; % 偏光板を0°～359°回転
        
        Ep = P * (Ba1*F_ice*Ba1') * [0; 1]; % 偏光状態(3.10)
        Es = P * (Ba1*F_ice*Ba1') * [1; 0]; % 偏光状態(3.11)

        T(t,1) = sum(percentage * abs(Es(:).^2) + (1-percentage) * abs(Ep(:).^2));
        
        X = X + (S(t,1) * T(t,1) * x1(t,1)); % 分子の積分計算
        Y = Y + (S(t,1) * T(t,1) * y1(t,1)); % 分子の積分計算
        Z = Z + (S(t,1) * T(t,1) * z1(t,1)); % 分子の積分計算
        
        kd = kd + (S(t, 1) * y1(t, 1)); % 分母の積分の計算 100は後で消えるので無視

    end

    X = X / kd;
    Y = Y / kd;
    Z = Z / kd;
    
    %---------XYZ座標をRGB変換----------
    % CIE XYZ系からCIE RGB系への変換
    R =     2.3655*X +  (-0.8971)*Y +  (-0.4683)*Z;
    G = (-0.5151)*X +      1.4264*Y +      0.0887*Z;
    B =     0.0052*X +  (-0.0144)*Y +      1.0089*Z;
    
    R = max(R, 0.0);
    G = max(G, 0.0);
    B = max(B, 0.0);

    %---------干渉色----------
    result1(:, ri, 1) = R;
    result1(:, ri, 2) = G;
    result1(:, ri, 3) = B;

end

%画像保存
imwrite(result1, ['./result_Rotation_ice5/iceThick--', num2str(ice_thick / 1000000), 'mm__degree--', num2str(angle_a), '.png']);
imshow(result1);
toc;