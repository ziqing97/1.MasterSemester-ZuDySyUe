function[y] = RK4_A5(y0,t0,a0,h,s)
l = (s-t0)/h;
t = zeros(l);
y = zeros(6,l);
t(1) = t0;
y(:,1) = y0;
a = a0;
for i=1:1:l
    k1 = dy6(t(i), y(:,i), a);
    k2 = dy6(t(i) + h/2, y(:,i) + k1 * h / 2, a);
    k3 = dy6(t(i) + h/2, y(:,i) + k2 * h / 2, a);
    k4 = dy6(t(i) + h, y(:,i) + h * k3, a);
    y(:,i+1) = y(:,i) + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
end