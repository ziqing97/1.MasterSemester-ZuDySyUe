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

wm = ones(2*n+1,1) / (2*n + 2*lambda);
wc = ones(2*n+1,1) / (2*n + 2*lambda);
wm(1) = lambda/(lambda+n);
wc(1) = wm(1) + (1 - alpha^2 + beta);

Chi1 = UT(x0,P0);

prae = praediktion(Chi1,dt);

xnnp = repmat(wm,1,2)'.* prae;
xnnp = sum(xnnp,2);

Pnnp = zeros(2);
for i = 1:2*n+1
    Pnnp = Pnnp + wc(i) * (prae(:,i)-xnnp) * (prae(:,i)-xnnp)';
end
Pnnp = Pnnp + Q;

Chi2 = UT(xnnp,Pnnp);

beob = beo(Chi2);
zn = beob;
zn = wm'.* zn;
zn = sum(zn);

S = 0;
for i = 1:2*n+1
    S = S + wc(i) * (beob(i)-zn) * (beob(i)-zn)';
end
R = 0.4;

S = S+R;

Pnnpp = zeros(2,1);
for i = 1:2*n+1
    Pnnpp = Pnnpp + wc(i) * (prae(:,i)-xnnp) * (beob(:,i)-zn)';
end


K = Pnnpp/S;

xnn = K*(15.4-zn) + xnnp;
Pnn = Pnnp - K*S*K';

function [xnnp] = praediktion(x,dt)
xnnp(1,:) = x(1,:) .* sin(x(1,:)) * dt;
xnnp(2,:) = x(2,:) .* cos(x(2,:)) * dt;
end

function [zn] = beo(x)
zn = x(1,:).*x(2,:) + x(1,:).^2;
end

function[Chi] = UT(x,P)
n = length(x);
alpha = 0.9;
kappa = 0;
Chi = zeros(n,length(x)*2+1);
Chi(:,1) = x;

lambda = alpha^2 * (n + kappa) - n;

Wurz = chol((n+lambda) * P);
Wurz = Wurz';

Chi(:,2:n+1) = repmat(x,1,2) + Wurz;
Chi(:,n+2:2*n+1) = repmat(x,1,2) - Wurz;
end