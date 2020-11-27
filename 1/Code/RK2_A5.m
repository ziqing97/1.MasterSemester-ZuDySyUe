function[y] = RK2_A5(y0,t0,a0,h,s)
l = (s-t0)/h;
t = zeros(l);
y = zeros(6,l);
t(1) = t0;
y(:,1) = y0;
a = a0;
for i=1:1:l
    k1 = dy6(t(i), y(:,i), a);
    k2 = dy6(t(i) + h, y(:,i) + k1 * h, a);
    y(:,i+1) = y(:,i) + h / 2 * (k1 + k2);
end
end