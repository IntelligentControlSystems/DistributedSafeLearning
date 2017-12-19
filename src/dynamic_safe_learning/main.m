% SAFE LEARNING
% Copyright (C) 2017 Andrea Carron
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

clear all
close all
clc

% DYNAMIC SAFE LEARNING FRAMEWORK: 
% From "Safe Learning for Distributed Systems with Bounded Uncertanties"
% The notation used in this software is the same of the paper


% initialize mpt for the plots
mpt_init;

% adding subfolders to the path
addpath(genpath('../common/functions'));
addpath(genpath('../common/scripts'));
addpath(genpath('./functions'));


% Systems definition
load_system;

% initial condition
x0 = [0.1; 0; 0; 0.1];

% tau variable for the optimization
tau     = [0.35;0.35];

% Dynamics Safe Sets computation
safe_sets = compute_dynamic_safe_sets(system, x0, tau);

% Simulation with Dynamics Safe Sets
x = x0;
for i = 1:100

    % applying random inputs
    u = 2*rand(2,1)-1;

    % checking if the input is safe
    old_input = [0;0];
    for j = 1:2
        % In this example the local state corresponds to the global
        safe = dynamics_safe_learning_input_check(system, safe_sets, x, u(j), j);
        if ~safe
            old_input(j) = u(j);
            u(j) = safe_sets.Li{j} * x;
        end
    end

    % generating the noise
    w = 2e-3 * rand(2,1) - 1e-3;
    while (w(1)'*system.Q(1,1)*w(1) > system.q(1) || ....
           w(2)'*system.Q(2,2)*w(2) > system.q(2) )
        w = 2e-3 * rand(2,1) - 1e-3;
    end

    % variables for plot
    y = sdpvar(2,1);
    E1 = YSet(y, y'*safe_sets.Pi{1}*y <= safe_sets.alphai{1} - x'*safe_sets.Gammai{1}*x);
    E2 = YSet(y, y'*safe_sets.Pi{2}*y <= safe_sets.alphai{2} - x'*safe_sets.Gammai{2}*x);
    future_state = system.A*x + system.B*u + system.G*w;

    % plot of system evolution
    figure(1)
    clf
    subplot(1,2,1)
    hold on
    E1.plot('wire',1)
    plot(x(1), x(2), 'kx', 'Markersize',10);
    plot(future_state(1),future_state(2), 'ko', 'Markersize',10);
    hold off
    axis square
    subplot(1,2,2)
    hold on
    E2.plot('wire',1)
    plot(x(3), x(4), 'kx', 'Markersize',10);
    plot(future_state(3),future_state(4), 'ko', 'Markersize',10);
    hold off
    axis square
    pause(0.1);

    % updating the system
    x = system.A*x + system.B*u + system.G*w;

end