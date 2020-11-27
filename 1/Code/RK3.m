function[y] = RK3(y0,t0,h,s)
l = (s-t0)/h;
t = zeros(1,l);
y = zeros(1,l);
t(1) = t0;
y(1) = y0;
for i=1:1:l
    k1=dy(t(i),y(i));
    k2=dy(t(i) + h/2, y(i) + k1 * h / 2);
    k3=dy(t(i) + h ,y(i) - h * k1 + 2 * h * k2);
    y(i+1) = y(i) + h / 6 * (k1 + 4 * k2 + k3);
    t(i+1) = t(i) + h;
end
end