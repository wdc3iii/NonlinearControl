function xdot = cartpole(x, u, params)
%PENDULUM Computes the continuous time dynamics of the pendulum
%   Inputs:
%   x:      state, (x, theta, xdot, thetadot). theta=0 corrosponds to upright
%   u:      input, (force on cart). positive force is right
%   params: struct, cart mass (M), pendulum mass (m), length (L), gravity (g)
% Outputs:
%   xdot: (xdot, thetadot, xddot, thetaddot).
xdot = cartpole_f(x, params) + cartpole_g(x, params) * u;
end