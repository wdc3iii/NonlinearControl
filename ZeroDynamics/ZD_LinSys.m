%% Setup (Linear system)
clear; clc; close all;

%% Easy System
A = [0 0; 1 1]; B = [1; 0];
Q = eye(2);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1, 1);
S12 = S(1, 2);

%% Plot Value Function
[y, z] = meshgrid(-5:0.1:5, -5:0.1:5);
v = y.^2 * S(1, 1) + 2 * y .* z * S(1, 2) + z.^2 * S(2, 2);
figure(1)
clf
surf(z, y, v)
xlabel('z')
ylabel('y')

%% Examine optimal trajectory vs. min Value restricted to z.
y0 = - S12 * z / S11;

figure(2)
clf
hold on
contour(z, y, v, linspace(0, 1000))
plot(reshape(z, 1, []), reshape(y0, 1, []), 'k')
xlabel('z')
ylabel('y')

z0IC = 3;
IC = [- S(1, 2) * z0IC / S(1, 1); z0IC];
[~, x] = ode45(@(t, x) A * x + B * (-K * x), [0, 10], IC);
plot(x(:, 2), x(:, 1), 'r')

[~, x] = ode45(@(t, x) A * x + B * ([0, -(1 - S12 / S11) * S12 / S11] * x), [0, 10], IC);
plot(x(:, 2), x(:, 1), 'g')
ylim([-5, 5])
hold off


%% Stabilize the black subspace
figure(2)
hold on
Kp = 5;
for r = -5:5
    for c = -5:5
        IC2 = [r; c];
        [~, x] = ode45(@(t, x) A * x + B * (-S12/S11 * [1 1] * x - Kp * [1 S12/S11] * x), [0, 10], IC2);
        plot(x(:, 2), x(:, 1), 'c')
    end
end
ylim([-5, 5])
hold off

%% Now consider a 3D System
A = [0 1 0; 0 0 0; 1 0 1]; B = [0; 1; 0];
eig(A)
Q = eye(3);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1:2, 1:2);
S12 = S(1:2, 3);

figure(3)
clf
hold on
z = linspace(-5, 5);
y = - S11 \ S12 * z;
plot3(z', y(1, :), y(2, :), 'k')
xlabel('z')
ylabel('y1')
zlabel('y2')

z0 = 3;
y0 = - S11 \ S12 * z0;
[~, x] = ode45(@(t, x) A * x + B * (-K * x), [0, 10], [y0; z0]);
plot3(x(:, 3), x(:, 1), x(:, 2), 'r')

C = [eye(2) S11 \ S12];
Kp = 5; Kd = 4;
[~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], [y0; z0]);
plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
hold off

%% Plot a bunch
figure(3)
hold on
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            IC2 = [r; c; d];
            [~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
        end
    end
end
hold off
eig(A - B * ([1 0]*C*A*A + Kp*[1 0]*C + Kd*[1 0]*C*A))

%% Now consider a 3D System, ZD depend on y2
A = [0 1 0; 0 0 0; 0 1 1]; B = [0; 1; 0];
eig(A)
Q = eye(3);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1:2, 1:2);
S12 = S(1:2, 3);

figure(4)
clf
hold on
z = linspace(-5, 5);
y = - S11 \ S12 * z;
plot3(z', y(1, :), y(2, :), 'k')
xlabel('z')
ylabel('y1')
zlabel('y2')

z0 = 3;
y0 = - S11 \ S12 * z0;
[~, x] = ode45(@(t, x) A * x + B * (-K * x), [0, 10], [y0; z0]);
plot3(x(:, 3), x(:, 1), x(:, 2), 'r')

C = [eye(2) S11 \ S12];
Kp = 5; Kd = 4;
% [~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], [y0; z0]);
[~, x] = ode45(@(t, x) A * x + B * ([1 2-Kp 2 - 2 * Kp] * x), [0, 20], [y0; z0]);
plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
hold off

%% Plot a bunch
figure(4)
hold on
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            IC2 = [r; c; d];
            [~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
        end
    end
end
hold off
eig(A - B * ([1 0]*C*A*A + Kp*[1 0]*C + Kd*[1 0]*C*A))

%% Stabilize a subspace y = -Kp * z
figure(6)
clf
figure(5)
clf
xlabel('z')
ylabel('y')
zlabel('doty')
hold on
A = [0 1 0; 0 0 0; 1 1 1]; B = [0; 1; 0];

% now stabilize the subspace y = -Kp * z
Kz = 5; % Nice stability
% Kz = 0.1; % Unstable
% Kz = 1.1; % Stable (surface barely stable)
% Kz = 0.9; % Unstable (surface barely unstable)
C = [1 0 Kz];

disp((1 - A(3, 2) / (1 + Kz * A(3,2)) * Kz) * (A(3,3) - A(3,1) * Kz))

[x, y] = meshgrid(-5:5, -5:5);
z = -x / Kz;
surf(z, x, y,'edgecolor','none')
shading interp
alpha(0.8)
hold off;
Kp = 10; 
Kd = 6;
Az = A(3, :);
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            figure(5)
            hold on
            IC2 = [r; c; d];
            [t, x] = ode45(@(t, x) A * x + B * 1/(1 + Kz * Az * B) * (-Kz * Az * A * x - Kp * C * x - Kd * (Kz * Az + [0 1 0]) * x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
            hold off

            figure(6)
            subplot(2,1,1)
            hold on
            plot(t, x(:, 1) + Kz * x(:,3))
            hold off
            subplot(2,1,2)
            hold on
            plot(t, (Kz * Az + [0 1 0]) * x')
            hold off
        end
    end
end
subplot(2,1,1)
xlabel('t')
ylabel('error from subspace y = -Kp z')

subplot(2,1,2)
xlabel('t')
ylabel('error dot from subspace y = -Kp z')

%% Examine zero dynamics
figure(7)
clf
hold on
[x, y] = meshgrid(-5:5, -5:5);
z = -x / Kz;
surf(z, x, y, 'FaceColor', [0.4940 0.1840 0.5560], 'edgecolor','none', 'FaceAlpha', 0.7)

z0 = linspace(min(z, [], 'all'), max(z, [], 'all'));
x0 = -Kz * z0;
y0 = (1 + Kz * A(3, 2)) \ (-Kz * (A(3,3) - Kz * A(3,1)) * z0);
plot3(z0, x0, y0, 'k', 'LineWidth', 1)
xlabel('z')
ylabel('y1')
zlabel('y2')

% Plot on zero dynamics
N = 10;
IC = [x0(N); y0(N); z0(N)];
[t, x] = ode45(@(t, x) A * x + B * 1/(1 + Kz * Az * B) * (-Kz * Az * A * x - Kp * C * x - Kd * (Kz * Az + [0 1 0]) * x), [0, 10], IC);
% plot3(x(:, 3), x(:, 1), x(:, 2), 'r', 'LineWidth', 2)

% Now minimize the value function
Q = eye(3);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1:2, 1:2);
S12 = S(1:2, 3);
z1 = linspace(-5, 5);
yV = - S11 \ S12 * z1;

y2 = -5:5;
N = size(y2, 2);
fimplicit3(@(z, y, dy) y + [1 0] * (S11 \ S12) * z, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor','none','FaceAlpha',.7)
plot3(z1, yV(1, :), yV(2, :), 'r')

x1 = -[1 0] * (S11 \ S12) * z1;
y1 = (1 + Kz * A(3, 2)) \ (-Kz * (A(3,3) - Kz * A(3,1)) * z0);
plot3(z1, x1, y1, 'c', 'LineWidth', 1)

legend('stabilizing $y = y_d$ surface', 'stable ZD surface', '$\min V|_z$ proj $y$', '$\min V|_z$','interpreter','latex')
set(gca,'FontSize',17)
% e = x(:, 1) + Kz * x(:, 3);
% de = x(:, 2) + (Kz * A(3, :) * x')';
% figure(8)
% clf
% plot(t, [e de])

%% Root Locusish ideas
A = [0 1 0; 0 0 0; 1 10 1]; B = [0; 1; 0];
figure(7)
clf
K = linspace(0, 10);
plot(K, (1 - A(3, 2) ./ (1 + K .* A(3,2)) .* K) .* (A(3,3) - A(3,1) .* K))
xlabel('K')
ylabel('a s.t. dotz = a z')
grid on


%% Compton System
k = 1; b = 0.2; m = 1;
% Nominal Dynamics with x = [x1 x2 x1dot x2dot]
Anom = [0 0 1 0;
     0 0 0 1;
     -k/m k/m -b/m b/m;
     k/m -k/m b/m -b/m];
Bnom = [0; 0; 1; 0];

% Choose the output collocated with actuator: y = x1, dy = x1dot, z1 = x2,
% z2 = x2dot
DPhi = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
A = DPhi * Anom / DPhi;
B = DPhi * Bnom;

% Feedback linearize output with u = - A(2, :)x + v
Ap = A - B * A(2, :);

% Segment matrices
Ay = Ap(3:4, 1);
Adoty = Ap(3:4, 2);
Az = Ap(3:4, 3:4);

% Design yd to stabilize the z dynamics
Q = eye(2);
R = 1;
[Kz, S, p] = lqr(Az, Ay, Q, R);
G = 10;
Kz = [G 2*sqrt(G)];

% Check if this stabilizes the whole system state
Az_feedback = ((1 - Adoty / (1 + Kz * Adoty) * Kz) * (Az - Ay * Kz));
eig(Az_feedback)

% Now stabilize the surface yd = -K * z
C = [1 0 Kz];
Kp = 5;
Kd = 4;
k_x = @(x) 1 / (C * A * B) * (- C * A * A - Kp * C - Kd * C * A) * x; 

% And of course try it
fh1 = figure(8); clf
fh2 = figure(9); clf
fh3 = figure(10); clf
z1 = 1;
z2 = 2;
for y1 = -5:5
    for y2 = -5:5
        Phi0 = [y1; y2; z1; z2];
        [t, x] = ode45(@(t, x) A * x + B * (-A(2, :) * x + k_x(x)), [0, 10], Phi0);

        figure(fh1)
        subplot(1,2,1)
        hold on
        plot(x(:, 1), x(:, 2), 'b')
        subplot(1, 2, 2)
        hold on
        plot(x(:, 3), x(:, 4), 'b')

        figure(fh2)
        subplot(2,1,1)
        hold on
        plot(t, C * x')
        subplot(2,1,2)
        hold on
        plot(t, C * A * x')

        figure(fh3)
        for i = 1:4
            subplot(4,1,i)
            hold on
            plot(t, x(:, i))
        end
    end
end

figure(fh1)
subplot(1,2,1)
xlabel('y')
ylabel('dy')
hold off
subplot(1,2,2)
xlabel('z1')
ylabel('z2')
hold off

figure(fh2)
subplot(2, 1, 1)
xlabel('t')
ylabel('e')
hold off
subplot(2, 1, 2)
xlabel('t')
ylabel('de')
hold off
