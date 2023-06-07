function f_x = cartpole_f(x, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (x, theta, xdot, thetadot). theta=0 corrosponds to upright
%   params: struct, cart mass (M), pendulum mass (m), length (L), gravity (g)
% Outputs:
%   f_x:    drift dynamics, f(x)
f_x = [
    x(3); 
    x(4);
    (-params.m * params.L * sin(x(2, :)) .* x(4, :).^2 + params.m * params.g * sin(x(2, :)) .* cos(x(2, :))) ...
    ./ (params.M + params.m * (1 - cos(x(2, :)).^2));
    (-params.m * params.L * sin(x(2, :)) .* cos(x(2)) .* x(4, :).^2 + (params.M + params.m) * params.g .* sin(x(2, :))) ...
    ./ ((params.M + params.m * (1 - cos(x(2, :)).^2)) * params.L);
];
end