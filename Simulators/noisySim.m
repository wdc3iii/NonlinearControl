function [t,x] = noisySim(model, controller, tspan, x0, params)
%CONTINUOUSSIM Performs a simulation with controller run 'continuously'
% Inputs
%   model:      function, takes state, input, model params
%   controller: function, takes state x as input (feedback controller)
%   tspan:      [t0, tf], bounds for simulation
%   x0:         initial condition for the simulation
%   params:     parameters for the model and the noise
[t, x] = ode45(@(t, x) model(x, controller(x), params), tspan, x0);
end

