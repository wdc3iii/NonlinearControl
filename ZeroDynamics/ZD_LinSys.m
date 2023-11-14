%% Setup (Linear system)
clear; clc; close all;

% % x = [y1 y2 z1 z2]'
% A = [0 1 0 0;
%     0 0 0 0;
%     0 0 0 1;
%     1 0 1 0];
% B = [0; 1; 0; 0];
% disp(rank(ctrb(A, B)) == 4)
% eig(A)
% Q = eye(4);
% R = 1;
% 
% sys = ss(A, B,[],[]);
% 
% % optimal feedback -Kx
% % value function V = x'Sx
% [K, S, P] = lqr(sys, Q, R);
% 
% %% Eval ICs
% z0 = [1; 1];
% eta_d = - S(1:2, 1:2) \ S(1:2, 3:4) * z0;
% 
% deta_d = - S(1:2, 1:2) \ S(1:2, 3:4) * A(3:4, :) * [eta_d; z0];

%% Easier System
A = [0 0; 1 1]; B = [1; 0];
Q = eye(2);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1, 1);
S12 = S(1, 2);

%% Plot Value Function
[y, z] = meshgrid(-5:0.1:5, -5:0.1:5);
v = y.^2 * S(1, 1) + 2 * y .* z * S(1, 2) + z.^2 * S(2, 2);
% figure()
% surf(z, y, v)
% xlabel('z')
% ylabel('y')

%% Examine optimal trajectory vs. min Value restricted to z.
y0 = - S12 * z / S11;

figure()
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


%% Stabilize the black subspace
Kp = 5;
for r = -5:5
    for c = -5:5
        IC2 = [r; c];
        [~, x] = ode45(@(t, x) A * x + B * (-S12/S11 * [1 1] * x - Kp * [1 S12/S11] * x), [0, 10], IC2);
        plot(x(:, 2), x(:, 1), 'c')
    end
end
ylim([-5, 5])

%% Now consider a 3D System
A = [0 1 0; 0 0 0; 1 0 1]; B = [0; 1; 0];
eig(A)
Q = eye(3);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1:2, 1:2);
S12 = S(1:2, 3);

figure()
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

%% Plot a bunch
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            IC2 = [r; c; d];
            [~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
        end
    end
end

eig(A - B * ([1 0]*C*A*A + Kp*[1 0]*C + Kd*[1 0]*C*A))

%% Now consider a 3D System, ZD depend on y2
A = [0 1 0; 0 0 0; 0 1 1]; B = [0; 1; 0];
eig(A)
Q = eye(3);
R = 1;

[K, S, P] = lqr(A, B, Q, R);
S11 = S(1:2, 1:2);
S12 = S(1:2, 3);

figure()
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

%% Plot a bunch
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            IC2 = [r; c; d];
            [~, x] = ode45(@(t, x) A * x + B * (-[1 0]*C*A*A*x - Kp*[1 0]*C*x - Kd*[1 0]*C*A*x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
        end
    end
end

eig(A - B * ([1 0]*C*A*A + Kp*[1 0]*C + Kd*[1 0]*C*A))

%% Stabilize a subspace y = -Kp * z
close all
figure(1)
xlabel('z')
ylabel('y')
zlabel('doty')
hold on
A = [0 1 0; 0 0 0; 1 1 1]; B = [0; 1; 0];

% now stabilize the subspace y = -Kp * z
% Kp = 5; % Nice stability
Kp = 0.1; % Unstable
% Kp = 1.1; % Stable (surface barely stable)
% Kp = 0.9; % Unstable (surface barely unstable)
C = [1 0 Kp];

[x, y] = meshgrid(-5:5, -5:5);
z = -x / Kp;
surf(z, x, y,'edgecolor','none')
shading interp
alpha(0.2)
hold off;
Kp1 = 10; 
Kd1 = 6;
Az = A(3, :);
for r = -5:2:5
    for c = -5:2:5
        for d = -5:2:5
            figure(1)
            hold on
            IC2 = [r; c; d];
            [t, x] = ode45(@(t, x) A * x + B * 1/(Kp * Az * B) * (-Kp * Az * A * x - Kp1 * C * x - Kd1 * (Kp * Az + [0 1 0]) * x), [0, 10], IC2);
            plot3(x(:, 3), x(:, 1), x(:, 2), 'b')
            hold off

            figure(2)
            subplot(2,1,1)
            hold on
            plot(t, x(:, 1) + Kp * x(:,3))
            hold off
            subplot(2,1,2)
            hold on
            plot(t, (Kp * Az + [0 1 0]) * x')

        end
    end
end
subplot(2,1,1)
xlabel('t')
ylabel('error from subspace y = -Kp z')

subplot(2,1,2)
xlabel('t')
ylabel('error dot from subspace y = -Kp z')


