%% Value Function Based Synthesis of Psi for Linear Systems
clear; clc; close all;

%% Define the system

% Define System
% Ay = [0; 1];
% Ady = [0; -2];
% Az = [0 1; -1 2];

% Compton System
load("ComptonSystem.mat")

% LQR
Q = diag([1 1 2 2]);
R = 0.1;

assert(size(Az, 1) == size(Az, 2))  % Az square (n-m x n-m)
assert(size(Az, 1) == size(Ay, 1))  % Ay has n-m rows
assert(size(Az, 1) == size(Ady, 1)) % Ady has n-m rows
assert(size(Ay, 2) == size(Ady, 2)) % Ay, Ady have m columns
m = size(Ay, 2);
n = m + size(Az, 1) / 2;
A = [zeros(m), eye(m), zeros(m, 2 * (n-m));
     zeros(m), zeros(m), zeros(m, 2 * (n-m));
     Ay, Ady, Az];
B = [zeros(m); eye(m); zeros(2*(n-m), m)];
if ~exist('Q', 'var')
    Q = eye(2 * n);
end
if ~exist('R', 'var')
    R = eye(m);
end

assert((size(Q, 1) == size(Q, 2)) && (size(Q, 1) == 2 * n)) % Q of proper size
assert((size(R, 1) == size(R, 2)) && (size(R, 1) == m)) % R of proper size

[Klqr, S, ~] = lqr(A, B, Q, R);
Sydy = S(1:m, m+1:2*m);
Syy = S(1:m, 1:m);
Syz = S(1:m, 2*m+1:end);
Sdydy = S(m+1:2*m, m+1:2*m);
Sdyz = S(m+1:2*m, 2*m+1:end);
Szz = S(2*m +1:end, 2*m+1:end);

Sbar = Sydy * (Sdydy - Sydy' * (Syy \ Sydy));
Psi = -Syy \ ((eye(m) + Sbar * Sydy' / Syy) * Syz - Sbar * Sdyz);

%% Question 1: Are the Zero Dynamics Stable?
assert(rank(eye(m) - Psi * Ady) == m)
Az_cl = Ay * Psi + Ady / (eye(m) - Psi * Ady) * Psi * (Ay * Psi + Az) + Az;
eig_Az_cl = eig(Az_cl);
if max(real(eig_Az_cl)) < 0
    disp("Zero Dynamics Stable")
else
    disp("Zero Dynamics Unstable")
end
disp(eig_Az_cl)

%% Design a feedback controller
% e = Cx
C = [eye(m) zeros(m), -Psi];
% u = -Kx
Kp = 10 * eye(m);
Kd = 2 * sqrt(10) * eye(m);
K = -(C * A * B) \ (C*A*A + Kp*C + Kd*C*A);

%% Value function still decreasing?
% Vdot = x'S(A + BK)x + (A + BK)'S
Vlqrdot = S * (A - B * Klqr) + (A - B * Klqr)' * S;
Vdot = S * (A + B * K) + (A + B * K)' * S;
[vlqr, dlqr] = eig(Vlqrdot);
[v, d] = eig(Vdot);

x0 = v(:, 4);

%% Simulate Response
% x0 = [1; 2; 3; 4];
% x0 = [-0.449518933872668;0.838072281062129;0.300382732593628;-0.073060206438843];
[t, x] = ode45(@(t, x) (A + B * K) * x, [0, 10], x0);

e = (C * x')';
de = ((C * A) * x')';

figure(1)
clf
subplot(3, 1, 1)
plot(t, [e de])
xlabel('$t$', 'interpreter', 'latex')
ylabel('$e$', 'interpreter', 'latex')
legend('$e$', '$\dot{e}$', 'interpreter', 'latex')
set(gca,'FontSize',17)

subplot(3, 1, 2)
plot(t, x(:, 3:4))
xlabel('$t$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
legend('$z_1$', '$z_2$', 'interpreter', 'latex')
set(gca,'FontSize',17)

val = diag(x * S * x');
valdot = diag(x * (S * (A + B * K) + (A + B * K)' * S) * x');

subplot(3, 1, 3)
plot(t, [val, valdot])
xlabel('$t$', 'interpreter', 'latex')
ylabel('$S$', 'interpreter', 'latex')
legend('$S(x)$', '$\dot{S}(x)$', 'interpreter', 'latex')
set(gca,'FontSize',17)




