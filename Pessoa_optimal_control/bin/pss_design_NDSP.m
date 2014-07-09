function pss_design_SP(system_name, verbose)

%
%   function pss_design_SP(system_name, N_sims, verbose)
%
%  
%   
%       INPUTS: system_name   - name of system for which this set will be
%                               used (used to inherit params_symb
%                               structure).
%               verbose       - Verbose level 0 (minimal), 1 (extra info), 2 (1+plots, just 2D systems), 3 (debug mode) 
%
%
%   Athanasios Tasoglou <tasoglou@gmail.com>, DCSC - TU Delft, 2013

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.4                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('------------------ Pessoa: Designing the NDSP Controller ---------------- ');


if verbose==3
	profile -memory on;
	setpref('profiler','showJitLines',1);
end


load(strcat(system_name, 'Controller_symb'),'params_symb');

pessoa_ndsp(system_name, verbose);


if verbose==3
	profile off;
	profsave(profile('info'),'profile_results/pss_build_cost');
end
