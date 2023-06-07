function [t,x,u] = continuousSim(model, controller, tspan, x0, params)
%CONTINUOUSSIM Performs a simulation with controller run 'continuously'
% Inputs
%   model:      function, takes state, input, model params
%   controller: function, takes state x as input (feedback controller)
%   tspan:   [t0, tf], bounds for simulation
%   x0:         initial condition for the simulation
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
if nargin(controller) == 1
    [t, x] = ode45(@(t, x) model(x, controller(x), params), tspan, x0, opts);
    u0 = controller(x(1, :)');
elseif nargin(controller) == 2
    [t, x] = ode45(@(t, x) model(x, controller(x, t), params), tspan, x0, opts);
    u0 = controller(x(1, :)', t(1));
end
u = zeros(length(u0), length(t));
for ii = 1:length(t)
    if nargin(controller) == 1
        u(:, ii) = controller(x(ii, :)');
    elseif nargin(controller) == 2
        u(:, ii) = controller(x(ii, :)', t(ii));
    end
end

