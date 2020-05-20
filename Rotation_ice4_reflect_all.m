%z�������ɕX�̎��������ƍl�������̃V�~�����[�V����
clear
tic; %���s���ԑ���

%--------------------------------�e�L�X�g�t�@�C���̓ǂݍ���----------------------------
%780nm����380nm�܂ł̐��l���L�������Ă����e�L�X�g�t�@�C��
x1=load('x-lambda.txt');
y1=load('y-lambda.txt');
z1=load('z-lambda.txt');
S=load('D65.txt');
%----------�o�͉摜��----------
height= 50;
width = 360;
result_all=zeros(height,width,3); % ���ʂ��\���摜��������  (�c,��,RGB)�̐}�̃T�C�Y
result2=zeros(height,width,3); 

%-------------�Ό��̐ݒ�---------------
rP = [0 pi/2 pi 0 0 0 0 0 0 0];% �Ό��̊p�x(10����
p=[cos(rP(1))^2 sin(rP(1))*cos(rP(1));sin(rP(1))*cos(rP(1)) sin(rP(1))^2];% �Ό��̍s���v�Z
p2=[cos(rP(2))^2 sin(rP(2))*cos(rP(2));sin(rP(2))*cos(rP(2)) sin(rP(2))^2];% �Ό��̍s���v�Z

%-------------�ʑ����t�B�����̐ݒ�-----------------
rB=[pi/4 pi*3/4 0 0 0 0 0 0 0 0];% �����ܐ��ޗ��̊p�x(10��)

Ba1=[cos(rB(1)) -sin(rB(1));sin(rB(1)) cos(rB(1))];% �����ܐ��ޗ��̉��]�v�Z
Ba2=[cos(rB(2)) -sin(rB(2));sin(rB(2)) cos(rB(2))];% �����ܐ��ޗ��̉��]�v�Z
Ba3=[cos(rB(3)) -sin(rB(3));sin(rB(3)) cos(rB(3))];% �����ܐ��ޗ��̉��]�v�ZR

%----------�ϑ�����------------------%
angle_a = 53;         %���ˊp(�u�����[�X�^�[�p:53.1��)
ice_thick = 500000; %�X�̌���[nm]
%���˕���(���ܗ�a->���ܗ�b)
%�X�̋��ܗ�(1.309,1.313)
n_air = 1;              %���C�̋��ܗ�(1)
n_x = 1.309;            %���˖ʂɐ��������̋��ܗ� x��
n_y = 1.313;            %y��
n_z = 1.309;            %z�� �����������ܗ�(c������)(1.313)
n_water = 1.333;        %���̋��ܗ�(1.333)

%----------�v�Z----------------------%
rad_a = deg2rad(angle_a);                   %���ˊp(���W�A��)
rad_b = asin( sin(rad_a) * (n_air/n_y) );     %�o�ˊp(���W�A��)
angle_b = rad2deg(rad_b);   %�o�ˊp(���܊p)(�u�����[�X�^�[�p:37.5��)
n = n_y * n_z / (sqrt(n_y.^2 * sin(rad_b).^2 + n_z.^2 * cos(rad_b).^2)); %���ˊp���ܗ�
a_dig = rad2deg(asin(sin(deg2rad(angle_a))/n_x));

rs_a2i = Rs(angle_a,n_air,n_x);
rp_a2i = Rp(angle_a,n_air,n_x);

rs_i2a = Rs(a_dig,n_x,n_air);
rp_i2a = Rp(a_dig,n_x,n_air);

rs_a2w = Rs(angle_a,n_air,n_water);
rp_a2w = Rp(angle_a,n_air,n_water);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% ここに湖底反射を加えてもあり  %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum_rps = rs_a2i+rp_a2i+rs_i2a+rp_i2a+rs_a2w+rp_a2w;

%percentage = rs / (rs + rp);   %s�Ό��̊���

%----------���^�� �� xyz �� xy,rgb ----------

ice = abs((n - n_x) * (ice_thick / cos(rad_b)));

% X�`Z �Ό����� �o�ˌ����x��������
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
    
    for t = 1:81% MATLAB�̓Y������1�����Ȃ̂ŁA���̂悤�ɂ��Ă���
        w = t*5+375;%�g��380�`780nm�܂�(������5[nm]�����ׂ�)
        
        F=[1 0;0 1]; 
        F_ice=[exp(-1i*pi*ice/w) 0;0 exp(1i*pi*ice/w)]; %�X�̔g�����Ƃ̍s���v�Z
        P=[cos(pi*(ri-1)/180)^2 sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180);sin(pi*(ri-1)/180)*cos(pi*(ri-1)/180) sin(pi*(ri-1)/180)^2];% �Ό��̍s���v�Z
        
        %a�`j�܂Ŏg���̂𓯂��悤�ɏ���
        %E=P1 *Ba1*F*Ba2 *Bb1*F4*Bb2 *Bc1*F*Bc2 *Bd1*F4*Bd2 *P2; %�o�ˌ����x�̌v�Z,�W���[���Y�}�g���b�N�X
        E_i = P * F * [1; 0];
        E2_i = P * F * [0; 1];
        E_a = P * (Ba1*F_ice*Ba1') * [1; 0];
        E2_a = P * (Ba1*F_ice*Ba1') * [0; 1];
        E_w = P * (Ba1*F_ice*Ba1') * [1; 0];
        E2_w = P * (Ba1*F_ice*Ba1') * [0; 1];
        
        I(t,1) =sum((rs_a2i/sum_rps) * abs(E_i(:).^2) +(rp_a2i/sum_rps) * abs(E2_i(:).^2)+(rs_i2a/sum_rps) * abs(E_a(:).^2) +(rp_i2a/sum_rps) * abs(E2_a(:).^2)+(rs_a2w/sum_rps) * abs(E_w(:).^2) +(rp_a2w/sum_rps) * abs(E2_w(:).^2));
      
        X=X+(S(t,1)*I(t,1)*x1(t,1));
        Y=Y+(S(t,1)*I(t,1)*y1(t,1)); % Y numerator ���q�̐ϕ��̌v�Z
        Z=Z+(S(t,1)*I(t,1)*z1(t,1)); % Z numerator ���q�̐ϕ��̌v�Z
        
        kd=kd+(S(t,1)*y1(t,1));% k denominator �����̐ϕ��̌v�Z 100�͌��ŏ������̂Ŗ���

    end

    X = X/kd;
    Y = Y/kd;
    Z = Z/kd;
    
    %kd�Ŋ����ĕۑ�
    xy(ri,1)=X/(X+Y+Z);% x���W
    xy(ri,2)=Y/(X+Y+Z);% y���W

    %----------XYZ���W��RGB�ϊ�----------
    % CIE XYZ�n����CIE RGB�n�ւ̕ϊ�
    R =    2.3655*X +  (-0.8971)*Y +  (-0.4683)*Z;
    G = (-0.5151)*X +     1.4264*Y +     0.0887*Z;
    B =    0.0052*X +  (-0.0144)*Y +     1.0089*Z;
    
    R = max(R,0.0);
    G = max(G,0.0);
    B = max(B,0.0);

    %----------���F----------
    result_all(:,ri,1)=R; % R�̒l
    result_all(:,ri,2)=G; % G�̒l
    result_all(:,ri,3)=B;% B�̒l
end

%----------�摜�ۑ�----------
imwrite(result_all, ['./result_reflect_all/',num2str(ice_thick/1000000),'mm_',num2str(angle_a),'.png']);
imshow(result_all);

toc; %���s���ԑ���
%**********END**********