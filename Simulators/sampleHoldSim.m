function [t,x,u] = sampleHoldSim(model, controller, tspan, x0, dt, params, holdLen)
%SAMPLEHOLDSIM Performs a simulation with controller run 'continuously'
% Inputs
%   model:      function, takes state, input, model params
%   controller: function, takes state x as input (feedback controller)
%   tspan:   [t0, tf], bounds for simulation
%   x0:         initial condition for the simulation
%   dt:         sample hold time
%   params:     parameters for the model
%   holdLen:    number of timesteps to evaluate sim during hold
if nargin <= 6
    holdLen = 10;
end
simSteps = ceil((tspan(2) - tspan(1)) / dt);
t = linspace(tspan(1), tspan(2), simSteps * holdLen + 1);
x = zeros(length(x0), simSteps * holdLen + 1);
u0 = controller(x0);
u = zeros(length(u0), simSteps * holdLen + 1);
x(:, 1) = x0;
for ii = 1:simSteps
    t0 = tspan(1) + (ii - 1) * dt;
    x0 = x(:, (ii-1)*holdLen + 1);
    u0 = controller(x0);
    [~, x_hold] = ode45(@(t, x) model(x, u0, params), linspace(t0, t0+dt, holdLen + 1), x0);
    x(:, (ii-1)*holdLen+1:ii*holdLen + 1) = x_hold';
    u(:, (ii-1)*holdLen+1:ii*holdLen + 1) = u0;
end
end

