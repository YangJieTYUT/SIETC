clear;close all;
% 读取序列图像
fileFolder = fullfile(pwd,'sequence','motion');
% fileFolder = fullfile(pwd,'sequence','6_3');
dirOutput = dir(fullfile(fileFolder,'motion*.512.tiff'));
% dirOutput = dir(fullfile(fileFolder,'6_3_*.tiff'));
fileNames = {dirOutput.name}';
F = numel(fileNames);

I = imread(fileNames{1});
[M,N] = size(I);
sequence = zeros([size(I) F],class(I));
sequence(:,:,1) = I;

for p = 2:F
    sequence(:,:,p) = double(imread(fileNames{p}));
end
[M1,N1,F] = size(sequence);
sequence = double(sequence);
threed_show(uint8(sequence))

% 测试彩色图像
% sequence = double(imread('lena_gray_512.tif'));
% [M1,N1,F] = size(sequence);

% 混沌密钥生成
[key1,key2] = key_gen(sequence,'SHA-512');


% 初始化测量矩阵
[X1,Y1] = EKPHM(key1,M1*N1);
CR = 0.5; t = 1; 
Phi = X1(1:CR * (M1/t) * (N1/t));
Phi = reshape(1 - 2 * Phi,CR * (M1/t), N1/t);
Phi = sqrt(2/(CR * M1/t)) * Phi;

% 生成Sbox
sbox = sbox_gen(key1,F);

% 对hash值进行一个gray变换以及ELGamal-Gray非对称加密。
% p=479;
% k_pr = 3;
% [public_key,k_E,y] = ELGenc(p,k_pr,hash_value);
% [re_x,re_hash_512] = ELGdec(p,k_E,k_pr,y);
% hash密钥非对称加解密完成

% 稀疏化
dwt_sequence = mydwt(sequence);
% threed_show(uint8(sequence))

% 置乱
[X2,Y2] = EKPHM([x0,y0,p0+alpha1,q0 + alpha2],2*M1*N1*F);
sc_seq = scramble(dwt_sequence,X2,Y2);
% threed_show(uint8(sc_seq))

% 压缩采样
y = zeros(CR*M1,N1,F);
for i = 1:F
    y(:,:,i) = stp(Phi,sc_seq(:,:,i));
end
% threed_show(uint8(y))

 % 量化
maxy=max(y(:));miny=min(y(:));
yq=(255*(y-miny))/(maxy-miny);
addition_key = mod(yq,1);
yq = floor(yq);
% threed_show(uint8(yq))

% S盒扩散
[M2,N2,F] = size(y);
[X3,Y3] = EKPHM(key2,M2*N2*F);
r = floor(rand(1)* 255);
C = sc_diffusion(yq,X3,Y3,sbox,r);

threed_show(uint8(C))

if M1== 256
    zlim([1,256])
    zticks([1,50,100,150,200,256])
    ylim([1,256])
    yticks([1,50,100,150,200,256])
end
if M1== 512
    zlim([1,512])
    zticks([1,100,200,300,400,512])
    ylim([1,512])
    yticks([1,100,200,300,400,512])
end