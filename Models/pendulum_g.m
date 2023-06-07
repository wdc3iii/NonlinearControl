function g_x = pendulum_g(x, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (theta, thetadot). theta=0 corrosponds to upright
%   params: struct, mass (m), length (L), gravity (g)
% Outputs:
%   g_x:    actuation matrix, g(x)
g_x = [
    0; 
    1 / (params.m * params.L^2)
];
end