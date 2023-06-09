function xdot = pendulum(x, u, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (theta, thetadot). theta=0 corrosponds to upright
%   u:      input, (torque). positive torque is CCW
%   params: struct, mass (m), length (L), gravity (g)
% Outputs:
%   xdot: (thetadot, thetaddot).
xdot = pendulum_f(x, params) + pendulum_g(x, params) * u;
end