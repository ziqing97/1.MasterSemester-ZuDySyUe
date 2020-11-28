%% Zustandsschätzung in dynamischen Systemen Übung 1
% Ziqing Yu 3218051
% 25/11/2020
clc
clear all 
close all

%% Aufgabe 1
% Erste Beobachtung y1 ist ein (40*1) Vektor
y1 = [2.68,2.96,3.66,4.66,1.43;
      2.73,2.97,3.52,4.81,1.82;
      2.63,3.11,3.16,4.48,1.52;
      3.10,2.91,3.47,4.56,1.58;
      2.78,3.06,3.41,4.66,1.72;
      2.79,3.07,3.42,4.61,1.62;
      2.78,3.11,3.36,4.58,1.56;
      2.80,3.07,3.47,4.56,1.62];
y1 = y1';
% Zweite Beobachtung y2
y2 = [2.81,3.12,3.33,4.61,1.56;
      2.78,3.06,3.38,4.56,1.61;
      2.86,3.03,3.39,4.57,1.72;
      2.80,3.76,4.22,5.31,1.62;
      2.81,3.87,4.08,5.36,1.56;
      2.78,3.81,4.13,5.31,1.61;
      2.86,3.78,4.14,5.32,1.72;
      2.80,3.76,4.22,5.31,1.62];
y2 = y2';

[x1,dx1,sigma1,e] = SeqAus(y1);
[x2,dx2,sigma2] = SeqAus(y2);


%% Aufgabe 3
s = 0.5;
h = 0.1;
y0 = 0;
t0 = 0;
% 3.Ordnung (Ergebnisse nach jeder Schrittweite werden dokumentiert)
y_3 = RK3(y0,t0,h,s);
% 4.Ordnung (Ergebnisse nach jeder Schrittweite werden dokumentiert)
y_4 = RK4(y0,t0,h,s);
% Differenz
diff34 = y_3 - y_4;

%% Aufgabe 4
% y_(n+m) = y_n + m * h * c, dann ist alles einfach

%% Aufgabe 5
% Anfangswert
yEP1 = [1.900041699219e7;1.696980371094e7;-1.411637695313e6;-3.696823120117e1;-2.619514465332e2;-3.568358421326e3];
aEP1 = [0;0;9.313225746155e-7];
tEP1 = 1000;

yEP2 = [1.826071386719e7;1.605853173828e7;-7.697145507813e6;-8.006420135498e2;-7.120389938354e2;-3.370163917542e3];
aEP2 = [-9.313225746155e-7;0;9.313225746155e-7];
tEP2 = 2800;

ts = 1900;

% Schrittweite
h1_vor = 100;
h1_rueck = -100;
h2_vor = 1;
h2_rueck = -1;
% 4.Ordnung
% Schrittweite 100s
y_100_vor = RK4_A5(yEP1,tEP1,aEP1,h1_vor,ts);
y_100_rueck = RK4_A5(yEP2,tEP2,aEP2,h1_rueck,ts);
y_100_vor_end = y_100_vor(:,end);
y_100_rueck_end = y_100_rueck(:,end);
% Schrittweite 1s
y_1_vor = RK4_A5(yEP1,tEP1,aEP1,h2_vor,ts);
y_1_rueck = RK4_A5(yEP2,tEP2,aEP2,h2_rueck,ts);
y_1_vor_end = y_1_vor(:,end);
y_1_rueck_end = y_1_rueck(:,end);

% 2.Ordnung
% Schrittweite 100s
y2_100_vor = RK2_A5(yEP1,tEP1,aEP1,h1_vor,ts);
y2_100_rueck = RK2_A5(yEP2,tEP2,aEP2,h1_rueck,ts);
y2_100_vor_end = y_100_vor(:,end);
y2_100_rueck_end = y_100_rueck(:,end);
% Schrittweite 1s
y2_1_vor = RK2_A5(yEP1,tEP1,aEP1,h2_vor,ts);
y2_1_rueck = RK2_A5(yEP2,tEP2,aEP2,h2_rueck,ts);
y2_1_vor_end = y_1_vor(:,end);
y2_1_rueck_end = y_1_rueck(:,end);