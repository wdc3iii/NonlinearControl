function xdot = pendulum(x, u, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (theta, thetadot). theta=0 corrosponds to upright
%   u:      input, (torque). positive torque is CCW
%   params: struct, mass (m), length (L), gravity (g)
% Outputs:
%   xdot: (thetadot, thetaddot).
xdot = [
    x(2); 
    -params.g / params.L * sin(x(1)) + u / (params.m * params.L^2)
];
end