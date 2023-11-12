%% Setup (Linear system)
clear;clc;close all;

% x = [y1 y2 z1 z2]'
A = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    1 0 1 0];
B = [0; 1; 0; 0];
disp(rank(ctrb(A, B)) == 4)
eig(A)
Q = eye(4);
R = 1;

sys = ss(A, B,[],[]);

% optimal feedback -Kx
% value function V = x'Sx
[K, S, P] = lqr(sys, Q, R);

[X, Y, U, W] = ndgrid(linspace(-1,1,10));

V = cell(size(X));
for i = 1:numel(X)
    % i
    [t, x] = ode45(@(t, x) A*x - B*K*x, [0 10],[X(i) Y(i) U(i) W(i)]);
    V{i} = x;
end

%% Plot all of the solutions projected onto the output space
% clf
% hold on
% for i = 1:numel(V)
%   plot(V{i}(:,3), V{i}(:,4));
%   drawnow
% end

% For linear systems, the value function is encoded by the feedback term K
% for nonlinear systems, we cannot construct sucha  feedback term, so we
% have to plan a single open loop trajectory, which when tracked leads to
% unstable behavior.
%% Value function

clf;

[t, x] = ode45(@(t, x) A*x + B*(-B'*(S*x)), [0 10],[1 1 2 2]);

plot(t,x)


%% Output feedback only is unstable (low level)

clf;
hold on

[t_d, x_d] = ode45(@(t, x) A*x - B*K*x, [0 10],[1 1 2 2]);

t_fine = interp1(1:length(t_d), t_d,linspace(1,length(t_d),100*length(t_d)));
x_d_fine = interp1(t_d, x_d, t_fine);

y = @(x, t) x(1) - x_d_fine(find(t_fine>=t,1),1);
dy = @(x, t) x(2) - x_d_fine(find(t_fine>=t,1),2);
ddy = @(t) [0 1 0 0]*(A-B*K)*x_d_fine(find(t_fine>=t,1),:)';

K_output = [10 2*sqrt(10)];

[t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(t)), [0 9],[1 1 2 2]);

plot(t_d, x_d)
hold on
plot(t,x)
%% Output feedback with replannig is stable (MPC)

clf;
hold on

K_output = [10 2*sqrt(10)];

timestep = 0.1;
t_end = 10;
T = 0;
X = [1 1 2 2];
while T(end) < t_end - timestep
    [t_d, x_d] = ode45(@(t, x) A*x - B*K*x, [T(end) T(end) + 10],X(end,:));

    t_fine = interp1(1:length(t_d), t_d,linspace(1,length(t_d),100*length(t_d)));
    x_d_fine = interp1(t_d, x_d, t_fine);

    y = @(x, t) x(1) - x_d_fine(find(t_fine>=t,1),1);
    dy = @(x, t) x(2) - x_d_fine(find(t_fine>=t,1),2);
    ddy = @(t) [0 1 0 0]*(A-B*K)*x_d_fine(find(t_fine>=t,1),:)';

    [t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(t)), [T(end) T(end)+timestep],X(end,:));
    T = [T; t(2:end)];
    X = [X; x(2:end,:)];
end

plot(T, X);

%% Centroidal MPC

clf;
hold on

A_z = A(3:4,3:4);
B_z = A(3:4,1:2);

Q = eye(2);
R = eye(2);

sys = ss(A_z, B_z,[],[]);

[K, S, P] = lqr(sys, Q, R);

K_output = [100 2*sqrt(100)];


timestep = 0.01;
t_end = 10;
T = 0;
X = [1 1 2 2];
X_D = [0 0];
while T(end) < t_end - timestep
    [t_d, x_d] = ode45(@(t, x) A_z*x - B_z*K*x, [T(end) T(end) + 10],X(end,3:4));

    n_d = (K*x_d')';
    n_d_1 = n_d(:,1);
    n_d_2 = gradient(n_d_1)./gradient(t_d);
    n_d_2_dot = gradient(n_d_2)./gradient(t_d);

    t_fine = interp1(1:length(t_d), t_d,linspace(1,length(t_d),100*length(t_d)));
    n_d_1_fine = interp1(t_d, n_d_1, t_fine);
    n_d_2_fine = interp1(t_d, n_d_2, t_fine);
    n_d_2_dot_fine = interp1(t_d, n_d_2_dot, t_fine);

    y = @(x, t) x(1) - n_d_1_fine(find(t_fine>=t,1));
    dy = @(x, t) x(2) - n_d_2_fine(find(t_fine>=t,1));
    ddy = @(t) n_d_2_dot_fine(find(t_fine>=t,1));

    [t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(t)), [T(end) T(end)+timestep],X(end,:));
    T = [T; t(2:end)];
    X = [X; x(2:end,:)];
    for j = 2:length(t)
        X_D = [X_D; [y(x(j,:), t(j)) dy(x(j,:), t(j))]];
    end

end

plot(T, X);

%% Output feedback with replannig on zero dynamics is sometimes(?) stable (ZD-MPC)


%%%%% WHAT IS HAPPENING
%%% WHY DOES THIS WORK AT ALL, AND WHY IS THERE SUCH A TIGHT BAND OF
%%% STABILITY?????????
clf;
hold on

% K_output = 0.8*[11 2*sqrt(10)];
Q = eye(4);
R = 1;

sys = ss(A, B,[],[]);

[K, S, P] = lqr(sys, Q, R);
K_output = K(1:2);

timestep = 0.1;
t_end = 10;
T = 0;
X = [1  1 2 2];
while T(end) < t_end - timestep
    [t_d, x_d] = ode45(@(t, x) A*x - B*K*x, [T(end) T(end) + 10],[0 0 X(end,3:4)]);

    t_fine = interp1(1:length(t_d), t_d,linspace(1,length(t_d),100*length(t_d)));
    x_d_fine = interp1(t_d, x_d, t_fine);

    y = @(x, t) x(1) - x_d_fine(find(t_fine>=t,1),1);
    dy = @(x, t) x(2) - x_d_fine(find(t_fine>=t,1),2);
    ddy = @(t) [0 1 0 0]*(A-B*K)*x_d_fine(find(t_fine>=t,1),:)';

    [t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(t)), [T(end) T(end)+timestep],X(end,:));
    T = [T; t];
    X = [X; x];
    % clf
    % plot(T, X)
    % hold on
    % plot(t_d, x_d)
    % pause
end

plot(T, X);

%% Map the zero dynamics to optimal output
% Fix the zero dynamics coordinate and plot the value function for values
% of the outputs

[Z1, Z2] = meshgrid(linspace(-10,10,100));
[X,Y] = meshgrid(linspace(-30,30,100));
OUTPUT = cell(size(Z1));
for k = 1:numel(Z1)
    Z = [Z1(k) Z2(k)];

    % minimizing the value function fixing z
    OUTPUT{k} = -S(1:2,1:2)\S(1:2,3:4)*[Z(1); Z(2)];

end

%% and then use for control via lookup (i.e. zero dynamics policies)
clf;

Q = eye(4);
R = 1;

sys = ss(A, B,[],[]);

[K, S, P] = lqr(sys, Q, R);
K_output = K(1:2);

timestep = 0.01;
t_end = 10;
T = 0;
X = [1  1 2 2];

OUTPUT_1 = cellfun(@(s) s(1), OUTPUT);
OUTPUT_2 = cellfun(@(s) s(2), OUTPUT);

X_D = [0 0];

while T(end) < t_end - timestep
    z = X(end,3:4);
    Y1 = interp2(Z1, Z2, OUTPUT_1, X(end,3), X(end,4));
    Y2 = interp2(Z1, Z2, OUTPUT_2, X(end,3), X(end,4));

    y = @(x, t) x(1) - [1 0]*(-S(1:2,1:2)\S(1:2,3:4)*x(3:4));
    dy = @(x, t) x(2) - [0 1]*(-S(1:2,1:2)\S(1:2,3:4)*x(3:4));
    ddy = @(x, t) [0 1 0 0]*(A-B*K)*[(-S(1:2,1:2)\S(1:2,3:4)*x(3:4))' X(end,3:4)]';

    [t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(x, t)), [T(end) T(end)+timestep],X(end,:));
    T = [T; t];
    X = [X; x];

    X_D = [X_D; repmat([Y1 Y2],length(t),1)];
end

T = T(2:end);
X = X(2:end,:);
X_D = X_D(2:end,:);

% subplot(2,2,[1 2])
hold on
plot(T, X);
plot(T, X_D)

% subplot(2,2,3)
% hold on
% plot3(X(:,3),X(:,4),X_D(:,1),'r','linewidth',6)
% surf(Z1, Z2, OUTPUT_1)
% subplot(2,2,4)
% hold on
% surf(Z1, Z2, OUTPUT_2)
% plot3(X(:,3),X(:,4),X_D(:,2),'r','linewidth',6)

%% or in MPC
clf;
hold on

% K_output = 0.8*[11 2*sqrt(10)];
Q = eye(4);
R = 1;

sys = ss(A, B,[],[]);

[K, S, P] = lqr(sys, Q, R);
K_output = K(1:2);

timestep = 0.1;
t_end = 10;
T = 0;
X = [1  1 2 2];

OUTPUT_1 = cellfun(@(s) s(1), OUTPUT);
OUTPUT_2 = cellfun(@(s) s(1), OUTPUT);

while T(end) < t_end - timestep
    z = X(end,3:4);
    Y_1 = interp2(Z1, Z2, OUTPUT_1, X(end,3), X(end,4));
    Y_2 = interp2(Z1, Z2, OUTPUT_2, X(end,3), X(end,4));
    
    [t_d, x_d] = ode45(@(t, x) A*x - B*K*x, [T(end) T(end) + 10],[Y_1 Y_2 X(end,3:4)]);

    t_fine = interp1(1:length(t_d), t_d,linspace(1,length(t_d),100*length(t_d)));
    x_d_fine = interp1(t_d, x_d, t_fine);

    y = @(x, t) x(1) - x_d_fine(find(t_fine>=t,1),1);
    dy = @(x, t) x(2) - x_d_fine(find(t_fine>=t,1),2);
    ddy = @(t) [0 1 0 0]*(A-B*K)*x_d_fine(find(t_fine>=t,1),:)';

    [t, x] = ode45(@(t, x) A*x + B*(-K_output*[y(x, t); dy(x, t)] + ddy(t)), [T(end) T(end)+timestep],X(end,:));
    T = [T; t];
    X = [X; x];
    % clf
    % plot(T, X)
    % hold on
    % plot(t_d, x_d)
    % pause
end

plot(T, X);