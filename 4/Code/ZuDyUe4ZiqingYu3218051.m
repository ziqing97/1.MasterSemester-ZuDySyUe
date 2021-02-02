%% Zustandsschaetzung in dynamischen Systemen Uebung 3
% Ziqing Yu 3218051
% 01/23/2021

clc
close all
clearvars

%% Aufgabe 1
sigma_r = sqrt(1);
sigma_p = 2;
P1 = zeros(1,201);
P1(1) = 100;

P2 = zeros(1,201);
P2(1) = 100;

F = 0;
G = 1;
A = [-F,G * eye(1) * sigma_p^2 * G';zeros(1),F'];
B = expm(A);
Phi = B(2,2)';
Q = Phi * B(1,2);
R = sigma_r^2;


for t = 1:200
    H = 0.5;
    Pnnp = Phi * P1(t) * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);
    P1(t + 1) = (1 - K * H) * Pnnp;
end

figure
plot(sqrt(P1))

for t = 1:200
    H = cos(1 + t/120);
    Pnnp = Phi * P2(t) * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);
    P2(t + 1) = (1 - K * H) * Pnnp;
end

hold on
plot(sqrt(P2))
legend('Instrument 1', 'Instrument 2')

%% Aufgabe 2
% vor
load('aufgabe2.mat');
beobachtung(:,3) = beobachtung(:,3) / 180 * pi; % deg to rad
P = 1;

sigma_zwd = sqrt(1e-8);
sigma_clk = sqrt(1e-12);
sigma_r = sqrt(0.001);

x1 = zeros(3,290);
x1(:,1) = [0.15;0;0];

F = [0,0,0;0,0,1;0,0,0];
G = [1,0,0;0,0,0;0,0,1];
W = diag([sigma_zwd^2,0,sigma_clk^2]);
A = [-F,G * W * G';zeros(3),F'] * 300;
B = expm(A);

Phi = B(4:6,4:6)';
Q = Phi * B(1:3,4:6);
R = sigma_r^2;

P3 = cell(290,1);
P3(1) = {eye(3) * P};

for t = 1:length(beobachtung)
    xnn_p = Phi * x1(:,t);
    H = [1/sin(beobachtung(t,3)), 1, 0];
    z = beobachtung(t,2);
    Pnnp = Phi * P3{t} * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);
    x1(:,t+1) = xnn_p + K * (z - H * xnn_p);
    P3(t+1) = {(eye(3) - K * H) * Pnnp};
end
% figure
% hold on
% plot(x1(1,:))
% plot(x1(2,:))
% title("vor")
% legend('zwd','clk')

% rueck
F = [0,0,0;0,0,1;0,0,0];
G = [1,0,0;0,0,0;0,0,1];
W = diag([sigma_zwd^2,0,sigma_clk^2]);
A = [-F,G * W * G';zeros(3),F'] * (-300);
B = expm(A);
x2 = zeros(3,290);
x2(:,1) = x1(:,end);
beobachtung2 = flip(beobachtung,1);
P4 = cell(290,1);
P4(1) = {eye(3) * P};
for t = 1:length(beobachtung)
    xnn_p = Phi * x2(:,t);
    H = [1/sin(beobachtung2(t,3)), 1, 0];
    z = beobachtung2(t,2);
    Pnnp = Phi * P4{t} * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);
    x2(:,t+1) = xnn_p + K * (z - H * xnn_p);
    P4(t+1) = {(eye(3) - K * H) * Pnnp};
end

x2 = flip(x2,2);

% figure
% hold on
% plot(x2(1,:))
% plot(x2(2,:))
% title('rueck')
% legend('zwd','clk')

% kombiniert
P4 = flip(P4);
P5 = cell(290,1);
x3 = zeros(3,290);
for t = 1:290
    P5(t) = {inv(inv(P3{t}) + inv(P4{t}))};
    x3(:,t) = P5{t} * (inv(P3{t}) * x1(:,t) + inv(P4{t}) * x2(:,t));
end
% figure
% hold on
% plot(x3(1,:))
% plot(x3(2,:))
% title("kombiniert")
% legend('zwd','clk')

figure
hold on
plot(x1(1,:))
plot(x2(1,:))
plot(x3(1,:))
legend('vor','rueck','komb')
title('zwd')

figure
hold on
plot(x1(2,:))
plot(x2(2,:))
plot(x3(2,:))
legend('vor','rueck','komb')
title('clk')

figure
hold on
plot(x1(3,:))
plot(x2(3,:))
plot(x3(3,:))
title('d-clk')
legend('vor','rueck','komb')

p1 = zeros(3,290);
p2 = zeros(3,290);
p3 = zeros(3,290);
for i = 1:290
    p1(1,i) = sqrt(P3{i}(1,1)); % zwd vor
    p1(2,i) = sqrt(P3{i}(2,2)); % clk vor
    p1(3,i) = sqrt(P3{i}(3,3)); % dclk vor
    
    p2(1,i) = sqrt(P4{i}(1,1)); % rueck
    p2(2,i) = sqrt(P4{i}(2,2));
    p2(3,i) = sqrt(P4{i}(3,3));
    
    p3(1,i) = sqrt(P5{i}(1,1)); % komb
    p3(2,i) = sqrt(P5{i}(2,2));
    p3(3,i) = sqrt(P5{i}(3,3));
end

figure
hold on
plot(p1(1,:))
plot(p2(1,:))
legend('vor','rueck')
title('std zwd')

figure
plot(p3(1,:))
title('std zwd kombiniert')

figure
hold on
plot(p1(2,:))
plot(p2(2,:))
legend('vor','rueck')
title('std clk')

figure
plot(p3(2,:))
title('std clk kombiniert')

figure
hold on
plot(p1(3,:))
plot(p2(3,:))
title('std d-clk')
legend('vor','rueck')

figure
plot(p3(3,:))
title('std d-clk kombiniert')