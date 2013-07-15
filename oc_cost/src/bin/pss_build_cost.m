function pss_build_cost(system_name, N_sims, verbose)

%
%   function pss_build_cost(system_name, N_sims, verbose)
%
%   Assigns a cost value to each state the system and the result is stored in an ADD. The cost of each state must be 
%	defined in the pss_cost_def.m. See online documentation on how to properly write the pss_cost_def.m script. 
%   
%       INPUTS: system_name   - name of system for which this set will be
%                               used (used to inherit params_symb
%                               structure).
%               N_sims        - size of batches to use in the construction
%               verbose       - Verbose level 0 (minimal), 1 (extra info), 2 (1+plots, just 2D systems), 3 (debug mode) 
%
%
%   Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>, DCSC - TU Delft, 2013

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.4                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('------------------- Pessoa: Creating Cost ADD Initiated ---------------- ');


if verbose==3
	profile -memory on;
	setpref('profiler','showJitLines',1);
end


filename_cost = strcat(system_name, 'Costs'); 

nbatch=N_sims;

load(strcat(system_name, '_symb'),'params_symb');

%
minc = params_symb.min(params_symb.xoind) - params_symb.min(params_symb.xoind); % or just zero :) 
maxc = params_symb.max(params_symb.xoind) - params_symb.min(params_symb.xoind);
%
totloops=ceil(prod((maxc-minc)+ones(size(maxc)))/nbatch);
%
disp('Creating the .add containing the Systems State Costs.');
pessoa_cost2add(filename_cost,minc,maxc,verbose);


% Create the Cost Adjacency Matrix and compute the APSP.
disp('Calculating the Deterministic APSP.');
pessoa_dsp(system_name, verbose);


if verbose==3
	profile off;
	profsave(profile('info'),'profile_results/pss_build_cost');
end
