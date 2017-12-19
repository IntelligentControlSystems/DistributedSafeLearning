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

function [Ai, Bi, Qi, Gi, qi] = global_2_local_system(system)
%GLOBAL_2_LOCAL_SYSTEM Extract all the local systems from the global one

% extracting the variables form the structure
[A,B,G,Q,q,~,~,~,~,n,m,p,~,~,neighbours] = system_2_variables(system);

% extracting the subsystems
cumulated_system_length = [0; cumsum(n)];
cumulated_input_length  = [0; cumsum(m)];
cumulated_noise_length  = [0; cumsum(p)];

for i=1:length(n)
    % creating the matrix to select the subsystem
    % left transformation (just identity in the correspondence of the agent
    % submatrix
    Txl = zeros(n(i), cumulated_system_length(end));
    Txl(:,cumulated_system_length(i)+1:cumulated_system_length(i+1)) = eye(n(i));
   
    % right transformation (identiy on the i-th row and the columns
    % representing a neigbours
    Txr = zeros(cumulated_system_length(end));
    for j = 1:length(neighbours{i})
        Txr(cumulated_system_length(j)+1:cumulated_system_length(j+1),...
           cumulated_system_length(j)+1:cumulated_system_length(j+1)) = eye(n(j));
    end
    
    % right matrix transformation for the input
    Tur = zeros(cumulated_input_length(end), m(i));
    Tur(cumulated_input_length(i)+1:cumulated_input_length(i+1),:) = eye(m(i));   
    
    % transformations
    Ai{i} = Txl * A * Txr;
    Bi{i} = Txl * B * Tur;
    Gi{i} = G(cumulated_system_length(i)+1:cumulated_system_length(i+1),...
              cumulated_noise_length(i)+1:cumulated_noise_length(i+1));
    Qi{i} = Q(cumulated_noise_length(i)+1:cumulated_noise_length(i+1),...
              cumulated_noise_length(i)+1:cumulated_noise_length(i+1));
    qi{i} = q(i);
end


end

