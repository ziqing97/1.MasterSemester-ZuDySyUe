% Dy6
clc
close all
clearvars

%
datafull = importdata('ue6_EKF.txt');
data = datafull.data;
data2 = flip(data);
clear datafull

% wahr
dt = 0.01;
t = 0.01:dt:12.5;
yw = 100 * cos(2*t-2);
xw = 100 * sin(t);
figure
plot(xw,yw)

% RW
% sigma_p = sqrt(0.1);
% sigma_r = sqrt(0.01);
% F = [0,0;0,0];
% G = [1,0;0,1];
% A = [-F,G * eye(2) * sigma_p^2 * G';zeros(2),F'] * dt;
% B = expm(A);
% Phi = B(3:4,3:4)';
% Q = Phi * B(1:2,1:2);
% R = sigma_r^2;
% 
% x1 = zeros(length(data)+1,2);
% x1(1,:) = [0,-41.61];
% x1 = x1';
% P1 = zeros(2,2,length(data)+1);
% P1(:,:,1) = [10,0;0,10];

xp1 = -130.56;
yp1 = 0;
xp2 = 13.10;
yp2 = -100.35;
% for i = 1:length(data)
%     xnn_p = Phi * x1(:,i);
%     para = sqrt((xp1 - x1(1,i))^2 + (yp1 - x1(2,i))^2);
%     z = data(i,2);
%     H = [(xp1 - x1(1,i))/para,(yp1 - x1(2,i))/para];
%     Pnnp = Phi * P1(:,:,i) * Phi' + Q;
%     K = Pnnp * H' * inv(H * Pnnp * H' + R);
%     x1(:,i+1) = xnn_p + K * (z - H * xnn_p);
%     P1(:,:,i+1) = (eye(2) - K * H) * Pnnp;
% end


% IRW
xa = 0;
ya = -41.61;
[x21,P21] = IRW(data,xa,ya,xp1,yp1);
[x22,P22] = IRW(data,xa,ya,xp2,yp2);



[x23,P23] = IRW(data2,x21(1,end),x21(3,end),xp1,yp1);
[x24,P24] = IRW(data2,x22(1,end),x22(3,end),xp2,yp2);

x23f = flip(x23,2);
P23f = flip(P23,3);

x24f = flip(x24,2);
P24f = flip(P24,3);

x25 = zeros(4,1251);
P25 = zeros(4,4,1251);
for i = 1:1251
    P25(:,:,i) = inv(inv(P21(:,:,i)) + inv(P23f(:,:,i)));
    x25(:,i) = P25(:,:,i) * (inv(P21(:,:,i)) * x21(:,i) + inv(P23f(:,:,i)) * x23f(:,i));
end

x26 = zeros(4,1251);
P26 = zeros(4,4,1251);
for i = 1:1251
    P26(:,:,i) = inv(inv(P22(:,:,i)) + inv(P24f(:,:,i)));
    x26(:,i) = P26(:,:,i) * (inv(P22(:,:,i)) * x22(:,i) + inv(P24f(:,:,i)) * x24f(:,i));
end

figure
hold on
plot(x23(1,:),x23(3,:))
plot(x24(1,:),x24(3,:))
plot(xw,yw)


figure
hold on
plot(x21(1,:),x21(3,:))
plot(x22(1,:),x22(3,:))
plot(xw,yw)

figure
hold on
plot(x25(1,:),x25(3,:))
plot(x26(1,:),x26(3,:))
plot(xw,yw)
