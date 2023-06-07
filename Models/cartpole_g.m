function g_x = cartpole_g(x, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (x, theta, xdot, thetadot). theta=0 corrosponds to upright
%   params: struct, cart mass (M), pendulum mass (m), length (L), gravity (g)
% Outputs:
%   xdot: (xdot, thetadot, xddot, thetaddot).
g_x = [
    0; 
    0;
    1 / (params.M + params.m * (1 - cos(x(2, :)).^2));
    cos(x(2, :)) / ((params.M + params.m * (1 - cos(x(2, :)).^2)) * params.L)
];
end