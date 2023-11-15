%% Simulation
clear; clc; close all;
%%% Parameters
l = 1;
g = -9.81;
mc = 1;
mp = 1;

%%% Double Pendulum Actuated in Middle (Acrobot)
syms x th dx dth real
q = [x; th];
dq = [dx; dth];
vars = [q; dq];
D = [mc+mp mp*l*cos(th); mp*l*cos(th) mp*l^2];
H = [-mp*l*dth^2*sin(th); mp*l*g*sin(th)];
B = [1;0];

f = [vars(3:4); -D \H];
g = [0; 0; D \ B];

N = [0 1];
y = B'*q;
dy = jacobian(y,vars)*f;
Lf2y = jacobian(dy,vars)*f;
LgLfy = jacobian(dy,vars)*g;
z1 = N*q;
z2 = N*D*dq;

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


syms n1_ n2_ z1_ z2_
z_vars = [n1_ n2_ z1_ z2_];

Phi_inv = [n1_; z1_; n2_; z2_ - n2_ * cos(z1_)]; % gives x from [n; z]
dynamics.Phi_inv = matlabFunction(Phi_inv,'vars',{z_vars});

assert(all(all(eval(jacobian(subs(Phi,vars,Phi_inv),z_vars)) == eye(4)))); % phi inverse is correct

omega = simplify(subs(omega_x, vars,Phi_inv));
f_hat = simplify(subs(f_hat_x, vars,Phi_inv));
g_hat = simplify(subs(g_hat_x, vars,Phi_inv));


Dzd = jacobian(omega,z_vars);
Lf_zd = Dzd * [f_hat; omega];
Lg_zd = Dzd * [g_hat; 0; 0];

dynamics.zd = matlabFunction(omega' ,'vars',{z_vars});
dynamics.Lf_zd = matlabFunction(Lf_zd ,'vars',{z_vars}); 
dynamics.Lg_zd = matlabFunction(Lg_zd ,'vars',{z_vars}); 

%% Simulate
dynamics.K_ll = [10 2*sqrt(10)];
dynamics.K_z = [10 -5];

tspan = [0, 20];

x0 = [1 1.5 1 0.4];
% x0 = [-7.5000   -0.5000    6.5587   -5.2558];

options = odeset('Events',@(t, x)explosionEvent(t, x));
[t,x] = ode45(@(t,x) CartpoleODE(t,x,dynamics), tspan, x0, options);

% Plot
figure(1);
clf
subplot(1, 2, 1)
plot(x(:, 1), x(:, 3))
xlabel('x')
ylabel('dx')
subplot(1, 2, 2)
plot(x(:, 2), x(:, 4))
xlabel('th')
ylabel('dth')

sgtitle('Phase Portraits')

figure(2);
clf
subplot(1, 2, 1)
plot(t, x(:, [1 3]))
xlabel('t')
ylabel('x, dx')
legend('x', 'dx')
subplot(1, 2, 2)
plot(t, x(:, [2 4]))
xlabel('t')
ylabel('th, dth')
legend('th', 'dth')
sgtitle('Time Series')

figure(3)
clf
e = dynamics.y(x')' + (dynamics.K_z * [dynamics.z1(x'); dynamics.z2(x')])';
de = dynamics.dy(x')' + (dynamics.K_z * dynamics.zd(dynamics.Phi(x')')')';

plot(t, [e, de])
legend('e', 'de')
xlabel('t')
ylabel('e, de')

%% Animation
figure(4)
clf
pause
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
e_fine = dynamics.y(x_fine')' + (dynamics.K_z * [dynamics.z1(x_fine'); dynamics.z2(x_fine')])';
de_fine = dynamics.dy(x_fine')' + (dynamics.K_z * dynamics.zd(dynamics.Phi(x_fine')')')';

plot(t_fine, [e_fine de_fine], 'LineWidth', 3)
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


%% Visualizing Zero Dynamics

dynamics.K_z = [10 10];

figure(6)
clf
hold on
N = 21;
zbar = zeros(2, N^2);
zdotbar = zeros(2, N^2);
z1s = linspace(-0.05, 0.05, N);
z2s = linspace(-0.05, 0.05, N);
max_e = 0;
for z1_ind = 1:N
    for z2_ind = 1:N
        disp(z2_ind + (z1_ind - 1) * N)
        % zdot = omega(y_d(z), dyd_dz * zdot, z)
        % Solve for zdot: 0 = zdot - omega(y_d(z), dyd_dz * zdot, z)
        zb = [z1s(z1_ind); z2s(z2_ind)];
        f = @(zdot) zdot - dynamics.zd([-dynamics.K_z * zb; -dynamics.K_z * zdot; zb]')';
        [zdot1, fv1, flag1] = fsolve(f, -zb);
        [zdot2, fv2, flag2] = fsolve(f, [0; 0]);
        [zdot3, fv3, flag3] = fsolve(f, -[0 1; -1 0] * zb);
        [zdot4, fv4, flag4] = fsolve(f, -[0 -1; 1 0] * zb);
        zbar(:, z1_ind + (z2_ind - 1) * N) = zb;
        zdotbar(:, z1_ind + (z2_ind - 1) * N) = zdot;

        if any([flag1, flag2, flag3, flag4] == 1)
            disp('yes')
        end


        % x0 = dynamics.Phi_inv([-dynamics.K_z * zb; -dynamics.K_z * zdot; zb]')';
        % [t,x] = ode45(@(t,x) CartpoleODE(t,x,dynamics), tspan, x0, options);
        % 
        % e = dynamics.y(x')' + (dynamics.K_z * [dynamics.z1(x'); dynamics.z2(x')])';
        % de = dynamics.dy(x')' + (dynamics.K_z * dynamics.zd(dynamics.Phi(x')')')';
        % max_e = max(max_e, max(abs(e)));
        % nz = dynamics.Phi(x')';
        % plot3(nz(:, 3), nz(:, 4), e, 'b')
    end
end

quiver(zbar(1, :), zbar(2, :), zdotbar(1, :), zdotbar(2, :), 'k')
xlabel('$z_1$', 'Interpreter', 'latex')
ylabel('$z_2$', 'Interpreter', 'latex')
hold off

fprintf('\nMaximum error from ZD surface: %e\n', max_e)


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
z_dot = d.zd(n_z')';
Lf_zd = d.Lf_zd(n_z');
Lg_zd = d.Lg_zd(n_z');

% Compute feedback linearizing input
u = (d.LgLfy(x) + K_z*Lg_zd) \ (-(d.Lf2y(x)+K_z*Lf_zd) - K_ll*[n(1) + K_z*z; n(2) + K_z*z_dot]);

% hack for plotting
u = max(min(u, 50), -50);
if abs(u) == 50
    % disp('..........Clipping..........')
end

% Apply input to the system
dx = d.f(x) + d.g(x)*u;
end

function [pos, isterminal,direction] = explosionEvent(t, x)
    pos = norm(x)-100;
    isterminal = 1;
    direction = 1;
end