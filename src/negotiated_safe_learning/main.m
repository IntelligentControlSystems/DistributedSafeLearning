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

% NEGOTIATED SAFE LEARNING FRAMEWORK: 
% From "Safe Learning for Distributed Systems with Bounded Uncertanties"
% The notation used in this software is the same of the paper


% initialize mpt for the plots
mpt_init;

% adding subfolders to the path
addpath(genpath('../common/functions'));
addpath(genpath('../common/scripts'));
addpath(genpath('../functions'));

% Systems definition
load_system;

% initial condition
x0 = [0; 0; 0; 0];

% optimization parameter
tau   = 0.6;

% Safe Sets computation
safe_sets = compute_negotiated_safe_sets(system, x0, tau);

% setting initial condition to state
x = x0;

for i = 1:100
    
    % generating random input
    u = randn(2,1);
    
    % finding if that input is feasible
    negotiated_output = negotiated_safe_learning_input_check(system, safe_sets, x, u);
    
    % CONCATENATE INPUT PROPERLY
    u_safe = [negotiated_output.input{1};negotiated_output.input{2}]    

    % generating the noise
    w = 2e-3 * rand(2,1) - 1e-3;
    while (w(1)'*system.Q(1,1)*w(1) > system.q(1) || ....
           w(2)'*system.Q(2,2)*w(2) > system.q(2) )
        w = 2e-3 * rand(2,1) - 1e-3;
    end

    % variables needed to generate the plot
    y = sdpvar(2,1);
    E1 = YSet(y, y'*safe_sets.Pi{1}*y <= negotiated_output.alphat{1});
    E2 = YSet(y, y'*safe_sets.Pi{2}*y <= negotiated_output.alphat{2});
    future_state = system.A*x + system.B*u_safe + system.G*w;
    future_state_no_safety = system.A*x + system.B*u + system.G*w;

    % plot of the system evolution
    figure(1)
    clf
    subplot(1,2,1)
    hold on
    E1.plot('wire',1)
    plot(x(1), x(2), 'kx', 'Markersize',10);
    plot(future_state(1),future_state(2), 'ko', 'Markersize',10);
    plot(future_state_no_safety(1),future_state_no_safety(2), 'ks', 'Markersize',10);
    hold off
    axis square
    subplot(1,2,2)
    hold on
    E2.plot('wire',1)
    plot(x(3), x(4), 'kx', 'Markersize',10);
    plot(future_state(3),future_state(4), 'ko', 'Markersize',10);
    plot(future_state_no_safety(3),future_state_no_safety(4), 'ks', 'Markersize',10);
    hold off
    legend('set','Current state', 'future safe state', 'future state without safety')
    axis square
    pause(2);

    % updating the system
    x = system.A*x + system.B*u_safe + system.G*w;
end
