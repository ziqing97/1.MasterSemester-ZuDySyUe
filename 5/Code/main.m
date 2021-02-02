% Dynamische System Uebung 5
% Ziqing Yu 3218051
% 01/02/2021

%% Initialisierung
clc
close all
clearvars

%% Aufgabe 2
x0 = [2;3];
P0 = [1,0.7;0.7,4];
dt = 1;
Q = [0.1,0.05;0.05,0.2];
z = 15.4;

alpha = 0.9;
beta = 2;
kappa = 0;

n = length(x0);
lambda = alpha^2 * (n + kappa) - n;

Chi = zeros(n,2*n+1);
Chi(:,1) = x0;
Chi(2:n+1) = x0 + chol((lambda + n) * P0);
Chi(2:n+1) = x0 - chol((lambda + n) * P0);