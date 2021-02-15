function [x2,P2] = IRW(data,xa,ya,xp1,yp1)
dt = 0.01;
sigma_p = sqrt(0.1);
sigma_r = sqrt(0.01);
F = zeros(4);
G = [0,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
A = [-F,G * eye(4) * sigma_p^2 * G';zeros(4),F'] * dt;
B = expm(A);
Phi = B(5:8,5:8)';
Q = Phi * B(1:4,1:4);
R = sigma_r^2;

x2 = zeros(length(data)+1,4);
x2(1,:) = [xa,0,ya,0];
x2 = x2';
P2 = zeros(4,4,length(data)+1);
P2(:,:,1) = diag([10,10,10,10]);
% Lineaisierung

for i = 1:length(data)
    xnn_p = Phi * x2(:,i);
    para = sqrt((xp1 - x2(1,i))^2 + (yp1 - x2(3,i))^2);
    z = data(i,2);
    H = [(xp1 - x2(1,i))/para,0,(yp1 - x2(3,i))/para,0];
    Pnnp = Phi * P2(:,:,i) * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);
    x2(:,i+1) = xnn_p + K * (z - H * xnn_p);
    P2(:,:,i+1) = (eye(4) - K * H) * Pnnp;
end
end