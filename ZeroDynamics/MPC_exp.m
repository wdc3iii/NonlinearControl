A = [0 1; 0 0];
B = [0; 1]; 
C = [0; 1];

Abar = [A B C; zeros(2, 4)];
T = 1.4;
x0 = [0; 0; -0.5; 1];

xf = [eye(2) zeros(2)] * expm(Abar * T) * x0;
[t, x] = ode45(@(t, x) A * x + B * x0(3) + C, [0, T], x0(1:2));
disp(xf - x(end, :)')