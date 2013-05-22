function pss_build_OOB_set(system_name,filename_set,verbose)

%
%   pss_build_OOB_set(params_symb,filename_set,verbose)
%
%   Constructs a set to be used to define a control problem for example. 
%   This function constructs a set consisting only of the Out-Of-Bounds state.
%   For the general function see "pss_build_set".
%   The set is saved in the form of a BDD.
%   
%       INPUTS: system_name - name of system for which this set will be
%                             used (used to inherit params_symb
%                             structure).
%               filename_set- name of file where the BDD is to be stored
%               verbose     - Verbose level 0 (minimal), 1 (extra info), 2 (1+plots, just 2D systems), 3 (debug mode) 
%
%
%   Manuel Mazo Jr. <m.mazo@rug.nl>, INCAS3-RuG 2010.

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.0                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('----------------------- Pessoa: Target Set Initiated ------------------- ');

load(strcat(system_name, '_symb'),'params_symb');

% transform target set to labels

% Target Set has to be Under-approximated
maxc=zeros(params_symb.n,1);
maxc(1)=params_symb.num(1)+1;
minc=maxc;

%%
totloops=1;
% Here we call the pss_target_set to find the target set
pessoa_set2bdd(filename_set,minc,maxc,verbose);

filename=strcat(filename_set,'_symb');
save(filename, 'params_symb');
  