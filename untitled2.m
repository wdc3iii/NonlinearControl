%% Stuff
clear; clc; close all;

%% Define systems

A = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    1 0 1 0];
B = [0; 1; 0; 0];

Q = eye(4);
R = 1;

sys = ss(A, B,[],[]);

[K, S, P] = lqr(sys, Q, R);

%% Examine zero dynamics

z0 = [1; 1];
yd = -0.5 * S(1:2, 1:2) \ S(1:2, 3:4) * z0;
dyd = -0.5 * S(1:2, 1:2) \ S(1:2, 3:4) * A(3:4, :) * [yd; z0];