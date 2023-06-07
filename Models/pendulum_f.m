function f_x = pendulum_f(x, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (theta, thetadot). theta=0 corrosponds to upright
%   params: struct, mass (m), length (L), gravity (g)
% Outputs:
%   f_x:    drift dynamics, f(x)
f_x = [
    x(2, :); 
    params.g / params.L * sin(x(1, :))
];
end