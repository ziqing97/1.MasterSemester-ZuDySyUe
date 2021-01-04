%% Zustandsschaetzung in dynamischen Systemen Uebung 3
% Ziqing Yu 3218051
% 21/12/2020

clc
close all
clearvars

%% a
% Wanderpunkten
Pyw = [5,1,4,2,4,3,5];
Pxw = [0.5,0.5,3,3,6,6,9];

% Wanderwinkeln
ww = zeros(length(Pxw) - 1,1);
wn = length(ww);
for i = 1 : wn
    ww(i) = atan2(Pyw(i+1) - Pyw(i), Pxw(i+1) - Pxw(i)); 
end


% Alle Punkten
% Abstand
s = 0.25;

% Linke Halbe
count =  1;
Pxl(1) = 0.5;
Pyl(1) = 5;
for i = 1 :wn
    while (sqrt((Pxl(count) - Pxw(i+1))^2 + (Pyl(count) - Pyw(i+1))^2) > 2*s) && (count < 5000)
        Pxl(count+1) = Pxl(count) + s * cos(ww(i));
        Pyl(count+1) = Pyl(count) + s * sin(ww(i));
        count = count + 1;
    end
end

% Rechte Halbe
Pxl = Pxl(1:end-2);
Pyl = Pyl(1:end-2);
Pxr = Pxl;
Pyl(end) = 5;
Pyr = 10 - Pyl;
Pxr = flip(Pxr);
Pyr = flip(Pyr);

% 
Px = [Pxl, Pxr(2:end)];
Py = [Pyl, Pyr(2:end)];

figure
scatter(Py,Px)
title('\sigma_r = 1')


sigma_r = 0.5;
% sigma_r = 1.0;
%% b: KF Random Walk
dt = 1;

F = [0,0;0,0];
sigma_p = 0.1;

G = [1,0;0,1];
A = [-F,G * eye(2) * sigma_p^2 * G';zeros(2),F'];
B = expm(A);
Phi = B(3:4,3:4)';
Q = Phi * B(1:2,3:4);
P = eye(2);


H = [1,0;0,1];

len = length(Px);
x_rw(:,1) = [5;0.5];
for i = 1:len
    xnn_p = Phi * x_rw(:,i);
    Pnn_p = Phi * P * Phi' + Q;
    z = [Py(i);Px(i)];
    R = eye(2) * sigma_r^2;
    K = Pnn_p * H' * inv(H * Pnn_p * H' + R);
    x_rw(:,i+1) = xnn_p + K * (z - H * xnn_p);
    P = (eye(2) - K * H) * Pnn_p;
end
hold on
plot(x_rw(1,:),x_rw(2,:),'k','Linewidth',2)



%% c: KF IRW
F = [0,1,0,0;0,0,0,0;0,0,0,1;0,0,0,0];
G = [0,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
A = [-F,G * eye(4) * sigma_p^2 * G';zeros(4),F'];
B = expm(A);
Phi = B(5:8,5:8)';
Q = Phi * B(1:4,5:8);
P = eye(4);
H = [1,0,0,0;0,0,1,0];
x_Irw(:,1) = [5;0;0.5;0];
for i = 1:len
    xnn_p = Phi * x_Irw(:,i);
    Pnn_p = Phi * P * Phi' + Q;
    z = [Py(i);Px(i)];
    R = eye(2) * sigma_r^2;
    K = Pnn_p * H' * inv(H * Pnn_p * H' + R);
    x_Irw(:,i+1) = xnn_p + K * (z - H * xnn_p);
    P = (eye(4) - K * H) * Pnn_p;
end
hold on
plot(x_Irw(1,:),x_Irw(3,:),'r','Linewidth',2)


for i = 1:131
   plot([x_rw(1,i+1),Py(i)],[x_rw(2,i+1),Px(i)],'c'); 
end
for i = 1:131
   plot([x_Irw(1,i+1),Py(i)],[x_Irw(3,i+1),Px(i)],'m'); 
end
legend('Beobachtung','RW','IRW')
%% Differenz
%
d_RWy = Py - x_rw(1,2:end);
d_RWx = Px - x_rw(2,2:end);
d_RW = sqrt(d_RWy.^2 + d_RWx.^2);

d_IRWy = Py - x_Irw(1,2:end);
d_IRWx = Px - x_Irw(3,2:end);
d_IRW = sqrt(d_IRWy.^2 + d_IRWx.^2);
figure
plot(d_RW)
hold on
plot(d_IRW)
title('Differenz')
legend('RW', 'IRW')