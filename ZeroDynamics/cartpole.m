%% Simulation
clear; clc; close all;

%%% Parameters
%%% Cartpole
l = 1;
gr = 9.81;
mc = 1;
mp = 1;
syms x th dx dth u real
q = [x; th];
dq = [dx; dth];
vars = [q; dq];
D = [mc+mp mp*l*cos(th); mp*l*cos(th) mp*l^2];
Dq = @(in)[mc+mp mp*l*cos(in(2, :)); mp*l*cos(in(2, :)) mp*l^2];
H = [-mp*l*dth^2*sin(th); -mp*l*gr*sin(th)];
B = [1; 0];

f = [vars(3:4); -D \H];
g = [0; 0; D \ B];

N = [0 1];


%%% Acrobot
% l=1;
% gr=9.81;
% m1=1;
% m2=1;
% lc = l/2;
% I1 = 1;
% I2 = 1;
% syms t1 t2 td1 td2 u real
% q = [t1; t2];
% dq = [td1; td2];
% vars = [q; dq];
% D = [I1+I2+m2*l^2+2*m2*l*lc*cos(t2) I2+m2*l*lc*cos(t2); I2+m2*l*lc*cos(t2) I2];
% Dq = @(in)[I1+I2+m2*l^2+2*m2*l*lc*cos(in(2,:)) I2+m2*l*lc*cos(in(2,:)); I2+m2*l*lc*cos(in(2,:)) I2];
% C = [-2*m2*l*lc*sin(t2)*td2 -m2*l*lc*sin(t2)*td2; m2*l*lc*sin(t2)*td1 0];
% G = -[m1*gr*lc*sin(t1)+m2*gr*(l*sin(t1)+lc*sin(t1+t2)); m2*gr*lc*sin(t1+t2)];
% H = C*[td1; td2] + G;
% B = [0;1];
%
% f = [vars(3:4); -inv(D)*H];
% g = [0;0; D \ B];
%
% N = [1 0];

y = B'*q;
dy = jacobian(y,vars)*f;
Lf2y = jacobian(dy,vars)*f;
LgLfy = jacobian(dy,vars)*g;
z1 = N*q;
z2 = N*D*dq;

%%% Dynamics
dynamics.f = matlabFunction(f,'vars',{vars});
dynamics.g = matlabFunction(g,'vars',{vars});
dynamics.y = matlabFunction(y,'vars',{vars});
dynamics.dy = matlabFunction(dy,'vars',{vars});
dynamics.Lf2y = matlabFunction(Lf2y,'vars',{vars});
dynamics.LgLfy = matlabFunction(LgLfy,'vars',{vars});
dynamics.z1 = matlabFunction(z1,'vars',{vars});
dynamics.z2 = matlabFunction(z2,'vars',{vars});
dynamics.n1 = dynamics.y;
dynamics.n2 = dynamics.dy;

Phi = [y; dy; z1; z2];
dPhi = jacobian(Phi, vars);

f_hat_x = [eye(2) zeros(2)]*dPhi*f;
g_hat_x = [eye(2) zeros(2)]*dPhi*g;
omega_x = [zeros(2) eye(2)]*dPhi*f;
assert(all(eval(simplify([zeros(2) eye(2)]*dPhi*g)) == [0; 0])); % zero dynamics had better not have control input.

dynamics.Phi = matlabFunction(Phi,'vars',{vars});
dynamics.dPhi = matlabFunction(dPhi,'vars',{vars});

syms n1_ n2_ z1_ z2_ real
z_vars = [n1_; n2_; z1_; z2_];

% Cartpole
% Phi_inv = [n1_; z1_; n2_; z2_ - n2_ * cos(z1_)]; % gives x from [n; z]

% Acrobot
% Phi_inv = [z1_; n1_; (z2_-(n2_*(cos(n1_)/2 + 1)))/(cos(n1_) + 3); n2_]; % gives x from [n; z]
Phi_inv = [B * n1_ + N' * z1_; B * n2_ + N' * ((N * Dq(B * n1_ + N' * z1_) * N') \ (z2_ - B' * Dq(B * n1_ + N' * z1_) * N' * n2_))];

dynamics.Phi_inv = matlabFunction(Phi_inv,'vars',{z_vars});

assert(all(all(eval(jacobian(subs(Phi,vars,Phi_inv),z_vars)) == eye(4)))); % phi inverse is correct

omega = simplify(subs(omega_x, vars, Phi_inv));
f_hat = simplify(subs(f_hat_x, vars, Phi_inv));
g_hat = simplify(subs(g_hat_x, vars, Phi_inv));


Domega = jacobian(omega, z_vars);
Lf_omega = Domega * [f_hat; omega];
Lg_omega = Domega * [g_hat; 0; 0];

dynamics.omega = matlabFunction(omega,'vars',{z_vars});
dynamics.Lf_omega = matlabFunction(Lf_omega ,'vars',{z_vars});
dynamics.Lg_omega = matlabFunction(Lg_omega ,'vars',{z_vars});
dynamics.K_ll = [10 2*sqrt(10)];
dynamics.K_z = [10 5];
options = odeset('Events',@(t, x)explosionEvent(t, x));

% For checking relative degree of error
LghatLfhate = simplify(subs(LgLfy, vars, Phi_inv) + dynamics.K_z*Lg_omega);
LgLfe = simplify(LgLfy + dynamics.K_z*subs(Lg_omega, z_vars, Phi));

% For MPC
df_dx = jacobian(f, vars);
dg_dxu = jacobian(g * u, vars);
dynamics.df_dx = matlabFunction(df_dx,'vars',{vars});
dynamics.dg_dxu = matlabFunction(dg_dxu,'vars',{vars, u});

%% Simulate
dynamics.K_ll = [10 2*sqrt(10)];
dynamics.K_z = [10 5];

tspan = [0, 20];

% x0 = [-0.200000000000000,0.320000000000000,2.253969504083596,-2.739547684553864];
x0 = [0.1, 0.1, 0, 0];

[t,x] = ode45(@(t,x) CartpoleODE(t,x,dynamics), tspan, x0, options);

%%% Animation
figure(4)
clf

set(gcf,'renderer','painters')
t_fine = linspace(t(1),t(end),1000);
x_fine = interp1(t,x,t_fine,'spline');

subplot(2,2,1)
plot(t_fine,[dynamics.y(x_fine'); dynamics.dy(x_fine')],'linewidth',3)
xlabel('$t$','interpreter','latex')
ylabel('$\eta$','interpreter','latex')
legend('$\eta_1$', '$\eta_2$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

subplot(2,2,2)
plot(t_fine,[dynamics.z1(x_fine'); dynamics.z2(x_fine')],'linewidth',3)
xlabel('$t$','interpreter','latex')
ylabel('$z$','interpreter','latex')
legend('$z_1$', '$z_2$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

subplot(2,2,4)
e_fine = dynamics.y(x_fine') + (dynamics.K_z * [dynamics.z1(x_fine'); dynamics.z2(x_fine')]);
de_fine = dynamics.dy(x_fine') + (dynamics.K_z * dynamics.omega(dynamics.Phi(x_fine')));

plot(t_fine, [e_fine; de_fine], 'LineWidth', 3)
xlabel('$t$','interpreter','latex')
ylabel('$e$','interpreter','latex')
legend('$e$', '$\dot{e}$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

subplot(2,2,3)
% axis off
base = rectangle('Position',[-0.2 -0.1 .4 0.2],'linewidth',2,'Curvature',0.2);
hold on;
P1_out = line([0 1],[0 1],'color','k','marker','o','markerSize',10,'linewidth',20);
P1_in = line([0 1],[0 1],'color',[0 0.4470 0.7410],'marker','o','markerSize',10,'linewidth',15);
axis([-10 10 -2 2])
axis equal
ax = gca;
set(ax,'YLim',[-1.5 1.5])
set(ax, 'XLim', [min(x(:, 1)) - 1, max(x(:, 1)) + 1])
buff = 5;
set(gca,'FontSize',17)
set(gca,'linewidth',2)

figure(18)
plot(t_fine, x_fine)

figure(4)
pause
tic
cont = true;
ind = 1;
while cont
    set(base,'Position',[-0.5+x_fine(ind,1) -0.375 1 .75])
    set(P1_out,'XData',[x_fine(ind,1) x_fine(ind,1)+sin(x_fine(ind,2))])
    set(P1_out,'YData',[0 cos(x_fine(ind,2))])
    set(P1_in,'XData',[x_fine(ind,1) x_fine(ind,1)+sin(x_fine(ind,2))])
    set(P1_in,'YData',[0 cos(x_fine(ind,2))])
    set(ax,'XLim',[x_fine(ind,1)-buff x_fine(ind,1)+buff])
    drawnow
    time = toc;
    ind = find(t_fine>time,1,'first');
    if isempty(ind)
        cont = false;
    end
end

figure()
plot(t_fine, x_fine)
%% Visualizing Zero Dynamics

figure(5)
clf
hold on
N = 21;
zbar = zeros(2, N^2);
zdotbar = zeros(2, N^2);
z1s = linspace(-0.4, 0.4, N);
z2s = linspace(-1, 1, N);
max_e = 0;
for z1_ind = 1:N
    for z2_ind = 1:N
        disp(z2_ind + (z1_ind - 1) * N)
        % zdot = omega(y_d(z), dyd_dz * zdot, z)
        % Solve for zdot: 0 = zdot - omega(y_d(z), dyd_dz * zdot, z)
        zb = [z1s(z1_ind); z2s(z2_ind)];
        % zb = z_fail;
        f = @(zdot) zdot - dynamics.omega([-dynamics.K_z * zb; -dynamics.K_z * zdot; zb]);
        [zdot1, fv1, flag1] = fsolve(f, -zb);
        [zdot2, fv2, flag2] = fsolve(f, [0; 0]);
        [zdot3, fv3, flag3] = fsolve(f, -[0 1; -1 0] * zb);
        [zdot4, fv4, flag4] = fsolve(f, -[0 -1; 1 0] * zb);
        zbar(:, z1_ind + (z2_ind - 1) * N) = zb;

        if any([flag1, flag2, flag3, flag4] == 1)
            if ~all([flag1, flag2, flag3, flag4] == 1)
                disp('here')
            end
            if flag1 == 1
                zdotbar(:, z1_ind + (z2_ind - 1) * N) = zdot1;
                zdot = zdot1;
                fv = fv1;
            elseif flag2 == 1
                zdotbar(:, z1_ind + (z2_ind - 1) * N) = zdot2;
                zdot = zdot2;
                fv = fv2;
            elseif flag3 == 1
                zdotbar(:, z1_ind + (z2_ind - 1) * N) = zdot3;
                zdot = zdot3;
                fv = fv3;
            else
                zdotbar(:, z1_ind + (z2_ind - 1) * N) = zdot4;
                zdot = zdot4;
                fv = fv4;
            end
            x0 = dynamics.Phi_inv([-dynamics.K_z * zb; -dynamics.K_z * zdot; zb]);
            [t,x] = ode45(@(t,x) CartpoleODE(t,x,dynamics), tspan, x0, options);

            e = dynamics.y(x') + (dynamics.K_z * [dynamics.z1(x'); dynamics.z2(x')]);
            de = dynamics.dy(x') + (dynamics.K_z * dynamics.omega(dynamics.Phi(x')));
            max_e = max(max_e, max(abs(e)));
            nz = dynamics.Phi(x')';

            plot3(nz(:, 3), nz(:, 4), e, 'b')
        else
            % Check if the output is valid relative degree here
            syms n2_temp real
            nz_temp = [-dynamics.K_z * zb; n2_temp; zb];
            x_temp = dynamics.Phi_inv(nz_temp);
            M = dynamics.LgLfy(x_temp) + dynamics.K_z * dynamics.Lg_omega(nz_temp);
        end
    end
end
quiver(zbar(1, :), zbar(2, :), zdotbar(1, :), zdotbar(2, :), 'k')
xlabel('$z_1$', 'Interpreter', 'latex')
ylabel('$z_2$', 'Interpreter', 'latex')
set(gca,'FontSize',17)

hold off

fprintf('\nMaximum error from ZD surface: %e\n', max_e)

%% Examine implicit zero dynamics, in a limiting sense:
N = 100;
z1_list = zeros(N, 1);
zdot_list = zeros(2, N);
z1_prev = 0;
z1_n = 0.2;
figure(6)
clf
hold on
zbar_final = [0;0];
zbardot_final = [0;0];
for z2_n = linspace(-0.5, 0.5, 31)
    for sgn = -1:2:1
        z1_prev = 0;
        z1_n = 0.2;
        for i = 1:N
            z1_list(i) = z1_n;
            zb = [z1_n; z2_n];
            f = @(zdot) zdot - dynamics.omega([-dynamics.K_z * zb; -dynamics.K_z * zdot; zb]);
            [zdot1, fv1, flag1] = fsolve(f, -zb);
            [zdot2, fv2, flag2] = fsolve(f, [0; 0]);
            [zdot3, fv3, flag3] = fsolve(f, -[0 1; -1 0] * zb);
            [zdot4, fv4, flag4] = fsolve(f, -[0 -1; 1 0] * zb);

            tmp = z1_n;
            if any([flag1, flag2, flag3, flag4] == 1)
                if flag1 == 1
                    zdot = zdot1;
                elseif flag2 == 1
                    zdot = zdot2;
                elseif flag3 == 1
                    zdot = zdot3;
                else
                    zdot = zdot4;
                end
                zdot_list(:, i) = zdot;
                plot(zb(1), zb(2), 'go')
                z1_n = z1_n + sgn * abs(z1_n - z1_prev);
                zbar_final = zb;
                zbardot_final = zdot;
            else
                plot(zb(1), zb(2), 'ro')
                z1_n = z1_n - sgn * abs(z1_n - z1_prev) / 2;
            end
            z1_prev = tmp;
        end
        x0 = dynamics.Phi_inv([-dynamics.K_z * zbar_final; -dynamics.K_z * zbardot_final; zbar_final]);
        [t,x] = ode45(@(t,x) CartpoleODE(t,x,dynamics), [0, 20], x0, options);
        nz = dynamics.Phi(x')';
        plot(nz(:, 3), nz(:, 4), 'b')
        M = dynamics.LgLfy(x0) + dynamics.K_z * dynamics.Lg_omega(dynamics.Phi(x0));
        plot3(zbar_final(1), zbar_final(2), M, 'ko')
        pause(0.01)
    end
end

%% Examine loss of relative degree
dynamics.LgLfe = matlabFunction(LgLfe ,'vars',{vars});
figure(7)
clf
fimplicit3(@(x,y,z) dynamics.LgLfe([0;x;y;z]))
xlabel('th')
ylabel('dx')
zlabel('dth')

dynamics.LghatLfhate = matlabFunction(LghatLfhate ,'vars',{z_vars});
figure(8)
clf
fimplicit3(@(x,y,z) dynamics.LghatLfhate([0;x;y;z]))
xlabel('n2')
ylabel('z1')
zlabel('z2')

%% Examine loss of relative degree
% dynamics.K_z = [3 0.8];
% Phi_bar = Phi - [-dynamics.K_z * [z1; z2]; -dynamics.K_z * dynamics.omega(dynamics.Phi(vars)')'; 0; 0];
% dPhi_bar = jacobian(Phi_bar, vars);
% charPoly = simplify(det(dPhi_bar));
% charPoly_f = matlabFunction(charPoly,'vars',{vars}');
%
% figure(7)
% fimplicit3(@(x,y,z) charPoly_f([0; x; y; z]), [-2*pi, 2*pi])
% xlabel('theta')
% ylabel('dx')
% zlabel('dth')

%% Zero dynamics surface
% Phi_bar = Phi - [-dynamics.K_z * [z1; z2]; -dynamics.K_z * dynamics.omega(dynamics.Phi(vars)')'; 0; 0];
% subs(Phi_bar,'x',0)
% functino = matlabFunction(Phi_bar(1),'vars',{vars}');
% figure(7)
% fimplicit3(@(x,y,z) functino([x; y; 0; z]), [-2*pi, 2*pi])
% xlabel('theta')
% ylabel('dx')
% zlabel('dth')

%% More
% x0 = [0; 0.4; 1; fzero(@(dth_) charPoly_f([0; 0.4; 1; dth_]), 0)];
%
% z_fail = [zeros(2) eye(2)] * dynamics.Phi(x0);

%% Quadratic shapes
% Define the positive definite matrix P
figure(22)
A = [0 1; -10, -2 * sqrt(10)];
P = lyap(A', eye(2));
% P = [2, 1; 1, 3];

% Define the range of x values
x1 = linspace(-0.2, 0.2, 100);
x2 = linspace(-0.6, 0.6, 100);

% Create a grid of x1 and x2 values
[X1, X2] = meshgrid(x1, x2);

% Initialize the Z values based on the quadratic form x' * P * x
Z = zeros(size(X1));

% Calculate the quadratic form for each point in the grid
for i = 1:numel(X1)
    x = [X1(i); X2(i)];
    Z(i) = x' * P * x;
end

% Plot the level curves
contour(X1, X2, Z, 20); % You can adjust the number of contour levels (e.g., 20)

% Add labels and title
xlabel('x1');
ylabel('x2');
title('Level Curves of x'' * P * x');

% Add a colorbar for reference
colorbar;

%% MPC Checking
T = 0.1;
x0 = zeros(4, 1);
u0 = 0;
f0 = dynamics.f(x0);
g0 = dynamics.g(x0);
df_dx0 = dynamics.df_dx(x0);
dg_dxu0 = dynamics.dg_dxu([x0; u0]);
Ac = df_dx0 + dg_dxu0;
Bc = g0;
Cc = f0 + g0 * u0 - Ac * x0 - Bc * u0;
Abar = [Ac Bc Cc; zeros(2, 6)];
eAT = expm(Abar * T);
Ad = eAT(1:4, 1:4);
Bd = eAT(1:4, 5);
Cd = eAT(1:4, end);

%% QP Checking
H = [2. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 2. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 1. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 1. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 2. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 2. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 1. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 1. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 1.];
f = [-0. -0. -0. -0. -0. -0. -0. -0.  0.];
A = [ 0.  0.  0.  0.  0.  0.  0.  0.  1.;
 -0. -0. -0. -0. -0. -0. -0. -0. -1.];
b = [ 50. 50.];
Aeq = [-1.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
0.00000000e+00;
  0.00000000e+00 -1.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   0.00000000e+00;
  0.00000000e+00  0.00000000e+00 -1.00000000e+00  0.00000000e+00 ...
   0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   0.00000000e+00;
  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.00000000e+00 ...
   0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   0.00000000e+00;
  1.00000000e+00 -4.98572308e-02  1.00000000e-01 -1.65111448e-03 ...
  -1.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 ...
   5.04114326e-03;
  0.00000000e+00  1.09971446e+00  0.00000000e+00  1.03302229e-01 ...
   0.00000000e+00 -1.00000000e+00  0.00000000e+00  0.00000000e+00 ...
  -5.08228652e-03;
  0.00000000e+00 -1.01339487e+00  1.00000000e+00 -4.98572308e-02 ...
   0.00000000e+00  0.00000000e+00 -1.00000000e+00  0.00000000e+00 ...
   1.01651114e-01;
  0.00000000e+00  2.02678973e+00  0.00000000e+00  1.09971446e+00 ...
   0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.00000000e+00 ...
  -1.03302229e-01];
beq = [-0.1 -0.1 -0.  -0.  -0.  -0.  -0.  -0. ];

sol = quadprog(H, f, A, b, Aeq, beq);

jaxOptSol = [ 9.99903934e-02  1.00131409e-01 -4.54981728e-05  4.87931482e-05 ...
  9.51522038e-02  1.09966991e-01 -9.84169912e-02  1.99844082e-01 ...
  3.07927186e-02];

jaxCVXPYSol = [ 1.00000000e-01  1.00000000e-01  4.49643191e-23 -1.38218301e-22 ...
  9.51692885e-02  1.09815169e-01 -9.82137866e-02  1.99502503e-01 ...
  3.07492943e-02];

osqpSol = [0.09998161646625116,0.0999865682484283,-2.0497113781387504e-06,5.443054055757252e-06,0.0951568010847993,0.10980617035751983,-0.09820623392361236,0.19949286881867795,0.030749307549233978];

cvxpySol = [ 1.00000000e-01  1.00000000e-01  4.49643191e-23 -1.38218301e-22,  9.51692885e-02  1.09815169e-01 -9.82137866e-02  1.99502503e-01,  3.07492943e-02];

disp(sol')
disp(jaxOptSol)
disp(jaxCVXPYSol)
disp(osqpSol)
disp(cvxpySol)

%% vis
sol = [ 1.00006191e-01,  1.12693039e-01,  1.43346055e-01...
         1.80902629e-01,  2.19980711e-01,  2.58136196e-01...
         2.94342948e-01,  3.28106937e-01,  3.58868279e-01...
         3.85458895e-01,  4.05409806e-01;
        1.00086626e-01,  9.23135520e-02,  7.55391037e-02...
         5.92631411e-02,  4.73066095e-02,  4.09524253e-02...
         4.06085763e-02,  4.67389583e-02,  6.05157790e-02...
         8.44793865e-02,  1.23489985e-01;
        9.58767306e-06,  2.54888891e-01,  3.60797302e-01...
         3.92879202e-01,  3.90527306e-01,  3.73520878e-01...
         3.50583324e-01,  3.23623124e-01,  2.89299861e-01...
         2.38566081e-01,  1.54073322e-01;
        2.25936145e-05, -1.59191382e-01, -1.82578521e-01...
        -1.48814214e-01, -9.46296799e-02, -3.48180983e-02...
         2.76176927e-02,  9.68586105e-02,  1.83088923e-01...
         3.03959199e-01,  4.88954527e-01];
u = [3.50511431 1.88403148 0.97906862 0.49463976 0.25785831 0.16550979 ...
  0.1531358  0.17579138 0.19399969 0.16008562];

%% More
e1 = [0.01; -0.5];
e1 = e1 / norm(e1);
l1 = -1/5;
e2 = [0.06; -0.2];
e2 = e2 / norm(e2);
l2 = -2;
A = [e1 e2] * [l1 0; 0 l2] / [e1 e2];
P = lyap(A', eye(2));

f = @(x,y) x.^2 * P(1,1) + 2 * x .* y * P(1,2) + y.^2 * P(2, 2);
[X,Y]=meshgrid(-0.2:0.001:0.2, -0.6:0.001:0.6);
z = f(X, Y);
figure(123)
clf
contour(X, Y, z, 100)
colorbar
[t, x] = ode45(@(t, x) A * x, [0, 20], [0.1; 0.3]);
hold on
plot(x(:, 1), x(:, 2), 'r')
%% Controller

function dx = CartpoleODE(t,x,d)
% Gains
K_ll = d.K_ll;
K_z = d.K_z;

% Grab eta, z coordinates
n_z = d.Phi(x);
n = [eye(2) zeros(2)] * n_z;
z = [zeros(2) eye(2)] * n_z;

% Compute z dynamics and derivatives
z_dot = d.omega(n_z);
Lf_omega = d.Lf_omega(n_z);
Lg_omega = d.Lg_omega(n_z);

% Compute feedback linearizing input
u = (d.LgLfy(x) + K_z*Lg_omega) \ (-(d.Lf2y(x)+K_z*Lf_omega) - K_ll*[n(1) + K_z*z; n(2) + K_z*z_dot]);

% hack for plotting
% u = max(min(u, 50), -50);
% if abs(u) == 50
%     % disp('..........Clipping..........')
% end

% Apply input to the system
dx = d.f(x) + d.g(x)*u;
end

function [pos, isterminal,direction] = explosionEvent(t, x)
pos = norm(x)-1000;
isterminal = 1;
direction = 1;
end