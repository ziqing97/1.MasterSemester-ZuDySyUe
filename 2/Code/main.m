%%
% Zustandsschaetzung in dynamischen Systemen Uebung 2
% Ziqign Yu 3218051

%% initial
clc
clearvars
close all

%% Aufgabe 1
syms dt

F1 = [0,0;0,0];
Phi1 = expm(F1 * dt);

F2 = [0,1,0,0;0,0,0,0;0,0,0,1;0,0,0,0];
Phi2 = expm(F2 * dt);

F3 = [0,1,0,0,0,0;
      0,0,1,0,0,0;
      0,0,0,0,0,0;
      0,0,0,0,1,0;
      0,0,0,0,0,1;
      0,0,0,0,0,0];
Phi3 = expm(F3 * dt);

%% Aufgabe 3
sigma = 1;
beta = 0.1;
F_A3 = -beta;
x = zeros(30,101);
x(:,1) = 0;
Phi_A3 = expm(F_A3 * 1);
for i = 1:30
    for t = 2:101
        x(i,t) = Phi_A3 * x(i,t-1) + sigma * randn(1);
    end
end
figure
hold on

varianz3 = zeros(1,30);
for i=1:30

    plot(x(i,:))
end
title('Aufgabe 3: 30 Realisation');
xlabel('Schritte')

% Varianz
for i = 1:100
    varianz3(i) = sqrt(var(x(:,i)));
end
figure
plot(varianz3.^2)
title('Varianz Aufgabe 3');
xlabel('Schritte')

%% Aufgabe 4
% Data Input
filename = 'E:\Studium\M1-Dynamischen Systemen\ZuDySy-Uebung\2\Code\random_a4.txt';
formatSpec = '%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
randoma4 = [dataArray{1:end-1}];
% clear inrelevant things
clearvars filename formatSpec fileID dataArray ans;

%
t = 1;
x_4 = zeros(100,21);
x_mittel = zeros(1,20);
x_4(:,1) = zeros(100,1);
F_A4 = [0,1;0,0];
Phi_A4 = expm(F_A4 * t);
for t = 2:21
    sum_mittel = 0;
    for j=1:50
        x_4((j-1)*2+1:2*j,t) = Phi_A4 * x_4((j-1)*2+1:2*j,t-1);
        x_4(2*j,t) = x_4(2*j,t) + randoma4(t-1,j);
        sum_mittel = sum_mittel + x_4((j-1)*2+1,t);
    end
    x_mittel(t-1) = sum_mittel/50;
end
figure
plot(x_mittel)
title('Aufgabe 4 Mittelwerte')
% Varianz

var4 = zeros(1,20);
for i=1:21
    var4(i) = var(x_4(:,i));
end
figure
plot(var4)

%% Aufgabe 5
omega = 0.1;
b = sqrt(2 * sqrt(2) * omega^3);
dt = 1;
F5 =  [0,1;-omega^ 2, -sqrt(2) * omega];
Phi5 = expm(F5 * dt);

G = [0;b];
Q = Phi5 * G * G' * Phi5';