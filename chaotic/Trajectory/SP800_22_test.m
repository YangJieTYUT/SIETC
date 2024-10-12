clc; clear; close all;

n = 1000000; 
% x0=0.6;y0=0.7;
% p=1;q=1;
pass_result = zeros(100,15);
result = zeros(100,15);

% for aaa = 1:100
% x0 = rand(1);y0 = rand(1);
% p=mod(rand(1) * 100,100);
% q=mod(rand(1) * 100,100);
x0 = 0.6;y0 = 0.9;
p=3.3;
q=5;
s = zeros(1,n);
t = 1000;h=0.001;
for i = 1:n  + t

    x1 = mod(exp(pi) * (((p * x0 * (1 - pi*y0) ^ 2)))+ q * y0 * x0 ^ 2 + pi * x0,1);
    y1 = mod(exp(pi) * (((p * y0 * (1 - pi*x0) ^ 2)))- q * x0 * y0 ^ 2 + pi * y0,1);
    x0 = x1;y0 = y1;
    if i > t
        %(x1+100)*2^16,mod 2
        s(i-t) = mod(floor((x1 + 100) * pow2(16)),2);
%                  s(i-t) = mod(floor((y1 + 100) * pow2(16)),2);
        %          s(i-t) = mod(floor((z1 + 100) * pow2(16)),2);
        %          s(i-t) = mod(floor((w1 + 100) * pow2(16)),2);
        %         if mod((i),3000) == 0
        if mod((i - t),3000) == 0
            %每隔3000次，加一次周期性小扰动
            x0 = x0 + h * sin(x0);
            %             y0 = y0 + h * sin(y0);
            %             z0 = z0 + h * sin(z0);
            %             w0 = w0 + h * sin(w0);
        end
    end
end


% 1.单比特测试  9725<k1<10275
clear K*; clear a b c h r t; clear x0 x1 y0 y1 z0 z1 w0 w1;
s1 = 2 * s - 1; sm = sum(s1);clear s1;sobs = abs(sm) / sqrt(n);

p_v01 = erfc(sobs / sqrt(2));
clear sm sobs;

% p_v01 = 0.7795;通过测试！！！

% 2.块内频率测试   每块序列为100，分成n/m块
m = 100; N = floor(n / m); f = zeros(1,N);
for i = 1:N
    f(i) = sum(s((i-1) * m + 1:i * m)) / m;
end

chai2 = 4 * m * sum((f - 1/2).^2);

p_v02 = gammainc(chai2/2,N/2,'upper');

% 感觉这个是正确的，正确的测试结果应该是0.7661
% p_v02 = gammainc(N/2,chai2/2,'upper');
clear m N f chai2;
%感觉这个是错误的
% 0.7661;0.2314;

% 3.游程测试

ps = sum(s)/n; r = bitxor(s(1:end - 1),s(2:end));vn_obs = 1 + sum(r);
p_v03 = erfc(abs(vn_obs - 2 * n * ps * (1 - ps)) / (2 * sqrt(2 * n) * ps * (1 - ps)));
clear ps r vn_obs;

% p_v03 = 0.4248;通过测试

%4.块内最长1游程测试

m = 10^4; K = 6; N = floor(n / m);

run1 = zeros(1,N);
% 分块计算
for i = 1:N
    s1 = s((i-1) * m + 1:i * m);
    j = 1; k = 1; t = 1;
    while j <= m
        if sum(s1(t:j)) == j - t + 1
            k = j - t + 1;j = j + 1;
        else
            t = t + 1; j = t + k;
        end
    end
    run1(i) = k;
end

%计算N个块中各类通常都的最长1游程的块数
v = zeros(1,7);
for i = 1:7
    if i ==1
        v(i) = sum(run1 < 11);
    elseif i == 7
        v(i) = sum(run1 > 15);
    else
        v(i) = sum(run1 == 9 + i);
    end
end
ps = [0.0882 0.2092 0.2483 0.1933 0.1208 0.0675 0.0727];
chai2_obs = sum((v - N * ps).^2./(N * ps));clear m N v run1 s1 i j k t ps;
p_v04 = gammainc(chai2_obs/2,K/2,'upper');
% p_v04 = gammainc(K/2,chai2_obs/2,'upper')

clear chai2_obs K;


% 5.二进制矩阵秩测试
M = 32;Q = 32; N =floor(100000/(M * Q));k0 = 0;k1 = 0;k2 = 0;
for i = 1:N
    R = zeros(M,Q);
    %找到N个M*Q的序列
    s1 = s((i-1) * M * Q + 1:i * M * Q);
    %   第j行的数字等于1到q，q+1到2q。。。。(j-1)q+1到jq
    for j = 1:M
        R(j,:) = s1((j-1) * Q + 1: j * Q);
    end
    %     初始化完成所有的矩阵
    
    %     开始求秩
    for a = 1:M - 1
        j = a + 1;
        while j <= M && R(a,a) == 0
            if R(j,a) == 1  %交换值
                t = R(j,:); R(j,:) = R(a,:); R(a,:) = t;
            end
            j = j + 1;
        end
        j = a + 1;
        while j <= Q && R(a,a) == 0
            if R(a,j) == 1
                t = R(:,j); R(:,j) = R(:,a); R(:,a) = t;
            end
            j = j + 1;
        end
        if R(a,a) == 1
            for u = a + 1:M
                R(u,:) = bitxor(R(u,:),R(a,:));
            end
            for u = a + 1:Q
                R(:,u) = bitxor(R(:,u),R(:,a));
            end
        end
    end
    r = sum(diag(R));
    
    if r == M
        k0 = k0 + 1;
    elseif r == M - 1
        k1 = k1 + 1;
    else
        k2 = k2 + 1;
    end
end

chai2_obs = (k0 - 0.2888 * N)^2/(0.2888 * N) + (k1 - 0.5776 * N)^2/(0.5776 * N) + ...
    (k2 - 0.1336 * N)^2/(0.1336 * N);
p_v05 = exp(-chai2_obs/2);
clear i j a chai2_obs k* M Q r R s1 t u N;

% 6.离散傅里叶测试(检测所有无周期)

n1 = 10^5; T = sqrt(log(1/0.05) * n1); N0 = 0.95 * n1 / 2;
X = 2 * s(1:n1) - 1;F = fft(X);F1 = abs(F(1:floor(n1/2)));
N1 = sum(F1 < T);d = (N1 - N0)/sqrt(0.95 * 0.05 * n1/4);
p_v06 = erfc(abs(d)/sqrt(2));
% p_v06 = 0.6845   通过测试
clear n1 T N0 N1 X F1 F d;

% 7.非重叠模板匹配测试

N = 20; M = floor(n/N); m = 10;
B = [1 0 0 1 1 0 1 0 1 1];%模板定长
u = (M - m + 1)/pow2(m);sigma2 = M * (1/pow2(m) - (2 * m - 1)/pow2(2 * m));
W = zeros(1,N);
for i = 1:N
    s1 = s((i-1) * M + 1:i * M);
    j = 1;
    while j <= M - m + 1
        if sum(bitxor(B,s1(j:j + m - 1))) == 0
            W(i) = W(i) + 1; j = j + m;
        else
            j = j + 1;
        end
    end
end

chai2_obs = sum((W - u) .^2/sigma2);
p_v07 = gammainc(chai2_obs/2,N/2,'upper');
% p_v07_1 = gammainc(N/2,chai2_obs/2,'upper')

clear B i j m M N chai2_obs u sigmafang s1 w;

% 8.重叠模板测试
N = 968; M = 1032; m = 9; B = [1 1 1 1 1 1 1 1 1];
lambta = (M - m + 1)/pow2(m); yita = lambta/2;
W = zeros(1,N);v = zeros(1,6);
for i = 1:N
    s1 = s((i - 1) * M + 1:i * M);
    j = 1;
    while j <= M - m + 1
        if sum(bitxor(B,s1(j:j + m - 1))) == 0
            W(i) = W(i) + 1;j = j + 1;
        else
            j = j + 1;
        end
    end
    if W(i)<5
        v(W(i) + 1) = v(W(i) + 1) + 1;
    else
        v(6) = v(6) + 1;
    end
end

ps = zeros(1,6);ps(1) = exp(-yita);
ps(2) = yita / 2 * exp(-yita);
ps(3) = yita * exp(-yita) / 8 * (yita + 2);
ps(4) = yita * exp(-yita) /8 * (yita .^ 2 / 6 + yita + 1);
ps(5) = yita * exp(-yita) /16 * (yita .^ 3/24 + yita .^ 2/2 + 3 * yita/2 + 1);
ps(6) = 1-sum(ps(1:5));
chai2_obs =sum((v - N * ps).^2./(N * ps));
% p_v08 = gammainc(5/2,chai2_obs/2,'upper')
p_v08 = gammainc(chai2_obs/2,5/2,'upper');
clear B i j m M N chai2_obs lambta yita s1 W ps v;


% 9.Maurer通用统计测试     检测序列能不能被显著压缩
% 两部分 Q*L; K*L;
L = 7; Q = 1280; K = floor(n/L) - Q;T = zeros(1,pow2(L));

for j = 1:Q
    t = s((j - 1) * L + 1:j * L);i = sum(t .* pow2(L - 1:-1:0));T(i + 1) = j;
end
sm = 0;
for j = Q + 1:Q + K
    t = s((j - 1) * L + 1:j * L);i = sum(t .* pow2(L - 1:-1:0));
    sm = sm + log2(j - T(i + 1));T(i + 1) = j;
end
c = 0.7 - 0.8/L + (4 + 32/L) * power(K, -3 * L)/15;
sigma = c * sqrt(3.125/K);
fn = sm/K;
p_v09 = erfc(abs((fn - 6.1962507)/(sqrt(2)*sigma)));

clear fn i j K L Q sm t T c;

% 10.线性复杂度测试   序列是否等价于使用长的LFSR（线性反馈移位寄存器）产生。
M = 1000; K = 6; N = floor(n/M);L = zeros(1,N);
for ii = 1:N
    S = s((ii - 1) * M + 1:ii * M);
    fx = zeros(1,N + 1);lx = zeros(1,N + 1);fx(1) = 1;lx(1) = 0;cond = 0;sm = 0;
    for i = 1:N
        d = mod(sum(fx(1:lx(i) + 1) .* S(i:-1:i-lx(i))),2);
        if d == 0
            lx(i + 1) = lx(i);
        else
            if sum(lx(1:i)) == 0
                lx(i + 1) = i;fx(lx(i + 1) + 1) = 1;
            else
                for j = i-1:-1:1
                    if lx(j) < lx(j + 1)
                        if j + 1 == i
                            cond = 1;break;
                        else
                            for k = i:-1:j + 1 + 1
                                if lx(k) == lx(k - 1)
                                    sm = sm + 1;
                                end
                            end
                            if sm == i - j - 1
                                sm = 0;cond = 1;break;
                            end
                        end
                    end
                end
                if cond == 1
                    cond = 0;
                    fx(i - j + 1:i - j + 1 + lx(j)) = mod(fx(i - j + 1:i - j + 1 + lx(j))...
                        +fx(1:1 + lx(j)), 2);
                    lx(i + 1) = max(lx(i),i - lx(i));
                end
            end
        end
    end
    L(ii) = max(lx);
end

u = M / 2 + (9 + (-1)^(M + 1))/36 - (M/3 + 2/9) / pow2(M);
T = (-1)^M * (L - u) + 2/9;
v = zeros(1,K+1);
for i = 1:N
    if T(i) <= -2.5
        v(1) = v(1) + 1;
    elseif T(i) <= -1.5
        v(2) = v(2) + 1;
    elseif T(i) <= -0.5
        v(3) = v(3) + 1;
    elseif T(i) <= 0.5
        v(4) = v(4) + 1;
    elseif T(i) <= 1.5
        v(5) = v(5) + 1;
    elseif T(i) <= 2.5
        v(6) = v(6) + 1;
    else
        v(7) = v(7) + 1;
    end
end
ps = [0.010417 0.03125 0.125 0.5 0.25 0.0625 0.020833];
chai2_obs = sum((v - N * ps) .^2./(N * ps));
p_v10 = gammainc(chai2_obs/2,K/2,'upper');
% p_v10_1 = gammainc(K/2,chai2_obs/2,'upper')


clear cond chai2_obs d fx i ii j k K L lx M u N S sigma sm T v ps;

% 11.序列测试   检验序列中长度为m的比特模板重复出现的次数

m = 3;
vm0 = zeros(1, pow2(m));vm1 = zeros(1,pow2(m-1));vm2 = zeros(1,pow2(m-2));
s1 = [s s(1:m-1)];
for i = 1:n
    for j = 1:pow2(m)
        bm0 = floor(mod((j - 1),2.^(m:-1:1))./(2.^(m-1:-1:0)));
        if sum(bitxor(bm0,s1(i:i + m - 1))) == 0
            vm0(j) = vm0(j) + 1;
        end
    end
    for j = 1:pow2(m-1)
        bm1 = floor(mod((j - 1),2.^(m-1:-1:1))./(2.^(m-2:-1:0)));
        if sum(bitxor(bm1,s1(i:i + m - 2))) == 0
            vm1(j) = vm1(j) + 1;
        end
    end
    for j = 1:pow2(m-2)
        bm1 = floor(mod((j - 1),2.^(m-2:-1:1))./(2.^(m-3:-1:0)));
        if sum(bitxor(bm1,s1(i:i + m - 3))) == 0
            vm2(j) = vm2(j) + 1;
        end
    end
end
psai2m0 = sum(vm0 .* vm0) * pow2(m)/n - n;
psai2m1 = sum(vm1 .* vm1) * pow2(m-1)/n - n;
psai2m2 = sum(vm2 .* vm2) * pow2(m-2)/n - n;

d_ps = psai2m0 -psai2m1;d_ps2 = psai2m0 - 2 * psai2m1 + psai2m2;

p_v11_1 = gammainc(d_ps,pow2(m-2),'upper');
p_v11_2 = gammainc(d_ps2,pow2(m-3),'upper');

% p_v11_1 = 0.1116;p_v11_2 = 0.0479;
clear b* d* i j m s1 v* ps*


% 12.近似熵测试
m = 3;
vm0 = zeros(1,pow2(m));vm1 = zeros(1,pow2(m + 1));
s1 = s(1:10000);n1 = length(s1);s1 = [s1 s1(1:m)];
for i = 1:n1
    for j = 1:pow2(m)
        bm0 = floor(mod((j - 1),2 .^(m:-1:1))./(2.^(m-1:-1:0)));
        if sum(bitxor(bm0,s1(i:i + m - 1))) == 0
            vm0(j) = vm0(j) + 1;
        end
    end
    for j = 1:pow2(m + 1)
        bm1 = floor(mod((j - 1),2 .^(m + 1:-1:1))./(2.^(m:-1:0)));
        if sum(bitxor(bm1,s1(i:i + m))) == 0
            vm1(j) = vm1(j) + 1;
        end
    end
end
vm0 = vm0/n1;vm1 = vm1/n1; fai0 = sum(vm0 .* log(vm0));fai1 = sum(vm1 .* log(vm1));
chai2 = 2 * n1 * (log(2) + fai1 - fai0);
p_v12 = gammainc(chai2/2,pow2(m - 1),'upper');
% p_v12_1 = gammainc(pow2(m - 1),chai2/2,'upper')
clear b* c* f* i j m n1 s1 v*;
% p_v12 = 0.4650;p_v12_1 = 0.4012;

% 13.累加和测试

n1 = 10000; s1 = s(1:n1); X = 2 * s1 - 1;
SL = zeros(1,n1);SR = zeros(1,n1);
SL(1) = X(1);SR(1) = X(end);
for i = 2:n1
    SL(i) = SL(i - 1) + X(1,i);SR(i) = SR(i - 1) + X(end - i + 1);
end
zL = max(SL); zR = max(SR);
kL1 = floor(4 * ( - n1/zL + 1)):floor(4 * (n1/zL - 1));
kL2 = floor(4 * ( - n1/zL - 3)):floor(4 * (n1/zL - 1));

p_v13L = 1 - sum(normcdf(zL * (4 * kL1 + 1)/sqrt(n1))...
    - normcdf(zL * (4 * kL1 - 1)/sqrt(n1)))...
    + sum(normcdf(zL * (4 * kL2 + 3)/sqrt(n1))...
    - normcdf(zL * (4 * kL2 + 1)/sqrt(n1)));
kR1 = floor(4 * ( - n1/zR + 1)); floor(4 * (n1/zR - 1));
kR2 = floor(4 * ( - n1/zR - 3)); floor(4 * (n1/zR - 1));

p_v13R = 1 - sum(normcdf(zR * (4 * kR1 + 1)/sqrt(n1))...
    - normcdf(zR * (4 * kR1 - 1)/sqrt(n1)))...
    + sum(normcdf(zR * (4 * kR2 + 3)/sqrt(n1))...
    - normcdf(zR * (4 * kR2 + 1)/sqrt(n1)));
clear i k* n1 s1 S* X z*;

% p_v13_L = 0.9992;p_v13_R = 1;


% 14.随机旅行测试

X = 2 * s - 1; S = zeros(1,n + 2);
for i = 2:n + 1
    S(i) = S(i - 1) + X(i - 1);
end
clear X;
v0 = zeros(1,8); v1 = v0;v2 = v0; v3 = v0;v4 = v0; v5 = v0;

j = 1;
for i = 2:n + 2
    if S(i) ~= 0
    else
        s1 = S(j:i); j = i;
        for x = -4:4
            if x == 0
            else
                if x < 0
                    id = x + 5;
                else
                    id = x + 4;
                end
                k = sum(s1 == x);
                if k == 0
                    v0(id) = v0(id) + 1;
                elseif k ==1
                    v1(id) = v1(id) + 1;
                elseif k ==2
                    v2(id) = v2(id) + 1;
                elseif k ==3
                    v3(id) = v3(id) + 1;
                elseif k ==4
                    v4(id) = v4(id) + 1;
                else
                    v5(id) = v5(id) + 1;
                end
            end
        end
    end
end

ps0 = zeros(1,8); ps1 = ps0; ps2 = ps0; ps3 = ps0; ps4 = ps0;
ps0(1:4) = 1 - 1./(2 * abs(-4:-1)); ps0(5:8) = fliplr(ps0(1:4));
ps1(1:4) = 1./(4 * (-4:-1).^2); ps1(5:8) = fliplr(ps1(1:4));
ps2(1:4) = 1./(4 * (-4:-1).^2).* (1 - 1./(2 * abs(-4:-1))); ps2(5:8) = fliplr(ps2(1:4));
ps3(1:4) = 1./(4 * (-4:-1).^2).* (1 - 1./(2 * abs(-4:-1))).^2; ps3(5:8) = fliplr(ps3(1:4));
ps4(1:4) = 1./(4 * (-4:-1).^2).* (1 - 1./(2 * abs(-4:-1))).^3; ps4(5:8) = fliplr(ps4(1:4));
ps5 = 1 - (ps0 + ps1 + ps2 + ps3 + ps4);
chai2_obs = zeros(1,8);
J = v0(1) + v1(1) + v2(1) + v3(1) + v4(1) + v5(1);
p_v14 = zeros(1,8);
for i = 1:8
    chai2_obs(i) = chai2_obs(i) + ...
        (v0(i) - J * ps0(i))^2/(J * ps0(i)) + ...
        (v1(i) - J * ps1(i))^2/(J * ps1(i)) + ...
        (v2(i) - J * ps2(i))^2/(J * ps2(i)) + ...
        (v3(i) - J * ps3(i))^2/(J * ps3(i)) + ...
        (v4(i) - J * ps4(i))^2/(J * ps4(i)) + ...
        (v5(i) - J * ps5(i))^2/(J * ps5(i));
    p_v14(i) = gammainc(chai2_obs(i)/2,5/2,'upper');
end
% p_v14
clear i id j J k ps* S s1 v* x chai2_obs
%     0.8955    0.8242    0.2926    0.4752
%     0.5731    0.1789    0.1094    0.0009  不通过！！！

% 15.随机变种测试
X = 2 * s - 1; S = zeros(1, n + 2);
for i = 2:n + 1
    S(i) = S(i - 1) + X(i - 1);
end
clear X;
J = sum(S == 0) - 1;
x = [-9:-1 1:9];
p_v15 = zeros(1,length(x));
for i=1:length(x)
    v=sum(S == x(i));
    p_v15(i) = erfc(abs(v-J)/(sqrt(2 * J * (4 * abs(x(i)) - 2))));
end
% p_v15
%
%     0.2283    0.4112    0.6687    0.9420
%     0.9104    0.8554    0.4504    0.3730
%     0.8850    0.0914    0.1477    0.0807
%     0.0361    0.0352    0.0338    0.0630
%     0.0904    0.0967    通过！！！！

clear i J S v x
% mean(p_v14)
% mean(p_v15)
p_value = [p_v01,p_v02,p_v03,p_v04,p_v05,p_v06,p_v07,p_v08,p_v09,p_v10,min(p_v11_1,p_v11_2),p_v12,min(p_v13L,p_v13R),min(p_v14),min(p_v15)];
% result(aaa,:) = p_value;
% pass_result(aaa,:) = p_value>0.01;
% end




