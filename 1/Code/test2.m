y1 = [2.68,2.96,3.66,4.66,1.43;
      2.73,2.97,3.52,4.81,1.82;
      2.63,3.11,3.16,4.48,1.52;
      3.10,2.91,3.47,4.56,1.58;
      2.78,3.06,3.41,4.66,1.72;
      2.79,3.07,3.42,4.61,1.62;
      2.78,3.11,3.36,4.58,1.56;
      2.80,3.07,3.47,4.56,1.62];
y = y1';

r = zeros(1,8);
for i = 1:8
    r(i) = 5 * i - 4; 
end
x = zeros(4,8);
dx = zeros(4,7);
sigma = zeros(1,8);
e = zeros(5,8);
% Design Matrix
A = [1,1,0,0;
     0,1,1,0;
     0,0,1,1;
     1,1,1,0;
     0,0,0,1];
% Gewicht Matrix
P = diag(ones(1,5) * (1/0.1^2));
x(:,1) = (A' * P * A) \ A' * P * y(:,1);
e(:,1) = y(:,1) - A * x(:,1);
sigma(1) = sqrt(e(:,1)' * P * e(:,1) / r(1)); 
Sigma = sigma(1)^2 * inv(A' * P * A);

for i = 2:8
    x(:,i) = x(:,i-1) + inv(sigma(i-1)^2 * inv(Sigma) + A' * P * A) * A' * P * (y(:,i) - A * x(:,i - 1));
    dx(:,i - 1) = x(:,i) - x(:, i-1);
    e(:,i) = y(:,i) - A * x(:,i);
    sigma(i) = sqrt((sigma(i-1)^2 * (r(i-1) + dx(:,i - 1)' * inv(Sigma) * dx(:,i - 1)) + e(:,i)' * P * e(:,i)) / r(i));
    Sigma = sigma(i)^2 * inv(sigma(i-1)^2 * inv(Sigma) + A' * P * A);
end