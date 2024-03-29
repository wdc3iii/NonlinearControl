%% Testing if Zero Dynamics Coordinates Preserved under transform
clear; clc; close all;

%% Collocated FBL
syms x th dx dth mc mp l g a1 a2 u
q = [x; th];        % robot state
dq = [dx; dth];     % robot velocity
% Dynamics components
D = [mc + mp, mp * l * cos(th); mp * l * cos(th), mp * l^2];
H = [-mp * l * dx * dth * sin(th); -mp * g * l * sin(th)];
B = [1; 0];
N = [0, 1];

% Construct coordinates
y = B' * q;
dy = B' * dq;
z1 = N * q;
z2 = N * D * dq;
% And construct the candidate diffeomorphism
I = [y; dy; z1; z2];
DI = jacobian(I, [q; dq]);
disp(DI)

%% Zero Dynamics Feedback
yd = @(a1, a2) -10 * a1 + 5 * a2;
% syms yd(a1, a2)
yf = y - yd(z1, z2);
dyf = simplify(dy - subs(jacobian(yd(a1, a2), [a1; a2]), {a1, a2}, {z1, z2}) * jacobian([z1; z2], [q; dq]) * [dq; D \ (B * u -H)]);

If = [yf; dyf; z1; z2];
DIf = simplify(jacobian(If, [q; dq]));
disp(DIf)
DIf = [DIf(1, :); DIf(3, :); DIf(2, :); DIf(4, :)];
DIf(1, :) = DIf(1, :) + subs(jacobian(yd(a1, a2), a2), {a1, a2}, {z1, z2}) * DIf(4, :);
DIf(1, :) = DIf(1, :) + subs(jacobian(yd(a1, a2), a1), {a1, a2}, {z1, z2}) * DIf(2, :);
disp(DIf)
det(subs(DIf(3:4, 3:4), {mc mp l g}, {1, 0.5, 0.5, 9.81}))
% solve(det(subs(DIf(3:4, 3:4), {mc mp l g}, {1, 0.5, 0.5, 9.81})))

%% Zero dynamics feedback, y_d (z1) only
% yd = @(a1) a1^3 + a1^2 + a1 + 3;
syms yd(a1)

yf = y - yd(z1);
dyf = simplify(dy - subs(jacobian(yd(a1), a1), a1, z1) * jacobian(z1, q) * dq);

If = [yf; dyf; z1; z2];
DIf = simplify(jacobian(If, [q; dq]));
disp(DIf)
DIf = [DIf(1, :); DIf(3, :); DIf(2, :); DIf(4, :)];
% DIf(1, :) = DIf(1, :) + subs(jacobian(yd(a1), a2), {a1, a2}, {z1, z2}) * DIf(4, :);
DIf(1, :) = DIf(1, :) + subs(jacobian(yd(a1), a1), {a1, a2}, {z1, z2}) * DIf(2, :);
disp(DIf)
det(subs(DIf(3:4, 3:4), {mc mp l g}, {1, 0.5, 0.5, 9.81}))

%% Evaluate full-rank-ness
% N = 50;
% thran = linspace(-pi, pi, N);
% 
% 
% max_cond = ones(N, 1);
% inds = ones(N, 3);
% parfor ii = 1:N
%     double(ii)
%     dxran = linspace(-10, 10, N);
%     dthran = linspace(-10, 10, N);
%     for jj = 1:N
%         for kk = 1:N
%             DIf_i = double(subs(DIf, {th dx dth mc mp l g}, ...
%                 {thran(ii), dxran(jj), dthran(kk), 1, 0.5, 0.5, 9.81}));
%             c = cond(DIf_i);
%             if rank(DIf_i) < 4
%                 disp('Rank deficient');
%             end
%             if c > max_cond(ii)
%                 max_cond(ii) = c;
%                 inds(ii, :) = [ii, jj, kk];
%             end
%         end
%     end
% end
% fprintf("Max Cond #: %e\n", max(max_cond));
