%% Compton System
clear; clc; close all;

k = 1; b = -2; m = 1;
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

Ay = Ap(3:4, 1);
Ady = Ap(3:4, 2);
Az = Ap(3:4, 3:4);

save('ComptonSystem', 'Ay', 'Ady', 'Az')