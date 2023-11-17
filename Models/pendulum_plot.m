function pendulum_plot(t,x,u, norm_bound, leg, formats)
%PENDULUM_PLOT Summary of this function goes here
%   Detailed explanation goes here
% And plotting
if nargin == 5
    [formats{1:length(t)}] = deal('-');
end
figure;
subplot(3, 1, 1)
hold on;
for ii = 1:size(t, 2)
    plot(t{ii}, x{ii}(1, :), formats{ii})
end
xlabel('Time (s)')
ylabel('$\theta$', 'Interpreter', 'latex')
hold off
subplot(3, 1, 2)
hold on;
for ii = 1:size(t, 2)
    plot(t{ii}, x{ii}(2, :), formats{ii})
end
xlabel('Time (s)')
ylabel('$\dot{\theta}$', 'Interpreter', 'latex')
hold off
subplot(3, 1, 3)
hold on;
for ii = 1:size(t, 2) - 1
    plot(t{ii}, vecnorm(x{ii} - x{end}), formats{ii})
end
plot(t{1}, norm_bound)
xlabel('Time (s)')
ylabel('$\|x - x_d\|$', 'Interpreter', 'latex')
legend(leg)
hold off
sgtitle('Pendulum State over Time')

figure;
hold on;
for ii = 1:size(t, 2) - 1
    plot(t{ii}, u{ii}, formats{ii})
end
hold off;
xlabel('Time (s)')
ylabel('Input (Nm)')
title('Pendulum Input over Time')

figure;
hold on;
for ii = 1:size(t, 2)
    plot(x{ii}(1, :), x{ii}(2, :), formats{ii})
end
hold off;
xlabel('$\theta$', 'Interpreter', 'latex')
ylabel('$\dot{\theta}$', 'Interpreter', 'latex')
legend(leg(1:end-1))
title('Pendulum Phase Plot')
end

