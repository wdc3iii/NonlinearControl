%% Attempt to debug the learning algorithm, for a linear system where I can compute everything by hand
clear; clc; close all;

% Setup
A = [0 1 0 0; 0 0 0 0; 0 0 0 1; 1 -2 -1 2];
B = [0 1 0 0];
Az = A(3:4, 3:4);
An1 = A(3:4, 1);
An2 = A(3:4, 2);

lam = 0.5;
P = [1.2514, 0.2771; 0.2771, 0.3166];
Klin = [0.8625 -0.4619];
nLin = [0.8335 1.9954];

syms z1 z2 th1 th2 th3 th4 real
z = [z1; z2];
th = [th1; th2; th3; th4];
psi1 = [th1 th2];
psi2 = [th3 th4];
psi = [psi1; psi2];

%% Compute stuff
w_z = simplify(Az * z + [An1 An2] * (psi * z));
residual_n2 = psi1 * w_z - psi2 * z;
V = simplify(z' * P * z);
Vdot = simplify(jacobian(V, z) * w_z);

Lv = Vdot + lam * V;
Lr = (residual_n2).^2;
Lz = max(0, Lv) + Lr;
gradL_th = jacobian(Lz, th);
v_grad = jacobian(Lv, th);
r_grad = jacobian(Lr, th);

func_grad = matlabFunction(gradL_th, 'vars', {[z', th']});
func_loss = matlabFunction(Lz, 'vars', {[z', th']});
func_loss_v = matlabFunction(Lv, 'vars', {[z', th']});
func_loss_r = matlabFunction(Lr, 'vars', {[z', th']});
func_vgrad = matlabFunction(v_grad, 'vars', {[z', th']});
func_rgrad = matlabFunction(r_grad, 'vars', {[z', th']});

%% Fix a current policy, th = bar_th. Examine gradient of loss as points z.
[z1_, z2_] = meshgrid(-3:0.1:3);
th1_ = Klin(1) + 0.5;
th2_ = Klin(2) - 0.1;
th3_ = 0;
th4_ = 0;
figure(1)
plot(Klin(1), Klin(2), 'ro')
th1s = th1_;
th2s = th2_;
th3s = th3_;
th4s = th4_;
lr = 0.005; 
for ii = 1:1000
    grad_1v = arrayfun(@(y1, y2) func_vgrad([y1, y2, th1_, th2_, th3_, th4_]) * [1; 0; 0; 0], z1_, z2_);
    grad_2v = arrayfun(@(y1, y2) func_vgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 1; 0; 0], z1_, z2_);
    grad_3v = arrayfun(@(y1, y2) func_vgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 0; 1; 0], z1_, z2_);
    grad_4v = arrayfun(@(y1, y2) func_vgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 0; 0; 1], z1_, z2_);
    grad_1r = arrayfun(@(y1, y2) func_rgrad([y1, y2, th1_, th2_, th3_, th4_]) * [1; 0; 0; 0], z1_, z2_);
    grad_2r = arrayfun(@(y1, y2) func_rgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 1; 0; 0], z1_, z2_);
    grad_3r = arrayfun(@(y1, y2) func_rgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 0; 1; 0], z1_, z2_);
    grad_4r = arrayfun(@(y1, y2) func_rgrad([y1, y2, th1_, th2_, th3_, th4_]) * [0; 0; 0; 1], z1_, z2_);

    loss = arrayfun(@(y1, y2) func_loss([y1, y2, th1_, th2_, th3_, th4_]), z1_, z2_);
    loss_v = arrayfun(@(y1, y2) max(0, func_loss_v([y1, y2, th1_, th2_, th3_, th4_])), z1_, z2_);
    loss_r = arrayfun(@(y1, y2) func_loss_r([y1, y2, th1_, th2_, th3_, th4_]), z1_, z2_);
    grad_1v(loss_v == 0) = 0;
    grad_2v(loss_v == 0) = 0;
    grad_3v(loss_v == 0) = 0;
    grad_4v(loss_v == 0) = 0;


    figure(2)
    surf(z1_, z2_, loss, 'EdgeColor','none')
    colorbar
    title('Loss')
    xlabel('z1')
    ylabel('z2')
    figure(3)
    surf(z1_, z2_, loss_v, 'EdgeColor','none')
    colorbar
    title('V Loss')
    xlabel('z1')
    ylabel('z2')
    figure(4)
    surf(z1_, z2_, loss_r, 'EdgeColor','none')
    colorbar
    title('R Loss')
    xlabel('z1')
    ylabel('z2')

    % figure(3)
    % quiver(z1_, z2_, grad_1, grad_2, 'k')
    % title('Gradient')
    % axis square
    % xlabel('z1')
    % ylabel('z2')
    % 
    th1_ = th1_ - lr * mean(grad_1v + grad_1r, 'all');
    th2_ = th2_ - lr * mean(grad_2v + grad_2r, 'all');
    th3_ = th3_ - lr * mean(grad_3v + grad_3r, 'all');
    th4_ = th4_ - lr * mean(grad_4v + grad_4r, 'all');
    % [~, inds] = max(loss, [], 'all');
    % th1_ = th1_ - lr * (grad_1v(inds) + grad_1r(inds));
    % th2_ = th2_ - lr * (grad_2v(inds) + grad_2r(inds));
    % th3_ = th3_ - lr * (grad_3v(inds) + grad_3r(inds));
    % th4_ = th4_ - lr * (grad_4v(inds) + grad_4r(inds));
    th1s = [th1s th1_];
    th2s = [th2s th2_];
    th3s = [th3s th3_];
    th4s = [th4s th4_];

    figure(1)
    clf
    subplot(1,2,1)
    hold on
    plot(Klin(1), Klin(2), 'ro')
    plot(th1s(1), th2s(1), 'go')
    plot(th1s, th2s, 'b.-')
    axis square
    legend('Linear Policy', 'Starting Point', 'Grad Descent')

    subplot(1,2,2)
    hold on
    plot(nLin(1), nLin(2), 'ro')
    plot(th3s(1), th4s(1), 'go')
    plot(th3s, th4s, 'b.-')
    axis square
    legend('Linear Policy', 'Starting Point', 'Grad Descent')
    hold off

    if max(loss, [], 'all') < 1e-10 && max(loss_v, [], 'all') == 0
        break
    end
    % pause(0.01)
end

%% Visualize loss at specific K
th1_ = 1;
th2_ = -100;
[z1_, z2_] = meshgrid(-3:0.005:3);
loss = arrayfun(@(y1, y2) func_loss([y1, y2, th1_, th2_]), z1_, z2_);
figure(7)
surf(z1_, z2_, loss, 'EdgeColor','none')
colorbar
title('Loss')
xlabel('z1')
ylabel('z2')

%% Visualize loss landscape (L2ish error)
[th1_, th2_] = meshgrid(-2:0.1:2);
% [th1_, th2_] = meshgrid(0.5:0.02:1, -0.5:0.02:0);
loss = arrayfun(@(t1, t2) mean(arrayfun(@(y1, y2) func_loss([y1, y2, t1, t2]), z1_, z2_), 'all'), th1_, th2_);
figure(5)
clf
hold on
surf(th1_, th2_, loss, 'EdgeColor', 'none')
% set(gca,'zscale','log')
plot3(Klin(1), Klin(2), 0.01, 'ro')
xlabel('theta_1')
ylabel('theta_2')
zlabel('loss')
hold off

%% Visualize loss landscape (Linfish error)
[th1_, th2_] = meshgrid(-2:0.1:2);
% [th1_, th2_] = meshgrid(0.5:0.02:1, -0.5:0.02:0);
loss = arrayfun(@(t1, t2) max(arrayfun(@(y1, y2) func_loss([y1, y2, t1, t2]), z1_, z2_),[], 'all'), th1_, th2_);
figure(6)
clf
hold on
surf(th1_, th2_, loss, 'EdgeColor', 'none')
% set(gca,'zscale','log')
plot3(Klin(1), Klin(2), 0.01, 'ro')
xlabel('theta_1')
ylabel('theta_2')
zlabel('loss')
hold off
