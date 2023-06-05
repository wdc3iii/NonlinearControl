function xdot = cartpole(x, u, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (x, theta, xdot, thetadot). theta=0 corrosponds to upright
%   u:      input, (force on cart). positive force is right
%   params: struct, cart mass (M), pendulum mass (m), length (L), gravity (g)
% Outputs:
%   xdot: (xdot, thetadot, xddot, thetaddot).
xddot = (-params.m * params.L * sin(x(2)) * x(4).^2 + params.m * params.g * sin(x(2)) * cos(x(2)) + u) ...
    / (params.M + params.m * (1 - cos(x(2)).^2));
thetaddot = (-params.m * params.L * sin(x(2)) * cos(x(2)) * x(4).^2 + (params.M + params.m) * params.g * sin(x(2)) + cos(x(2)) * u) ...
    / ((params.M + params.m * (1 - cos(x(2)).^2)) * params.L);
xdot = [
    x(3); 
    x(4);
    xddot;
    thetaddot
];
end