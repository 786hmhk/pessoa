function pss_build_set(system_name,target_set,filename_set, N_sims,type,verbose)

%
%   pss_build_set(params_symb,target_set,filename_set,N_sims,verbose)
%
%   Constructs a set to be used to define a control problem for example. The set is saved in the form of a BDD.
%   
%       INPUTS: system_name - name of system for which this set will be
%                             used (used to inherit params_symb
%                             structure).
%               target_set  - desired target set, when type=0.
%                             n-dim paralepiped containing the target set if
%                             type=1 (target_set=xset always works).
%               filename_set- name of file where the BDD is to be stored
%               N_sims      - size of batches to use in the construction
%               type        - 0: build n-dim paralepiped,
%                             1: pss_target_set.m
%               verbose     - Verbose level 0 (minimal), 1 (extra info), 2 (1+plots, just 2D systems), 3 (debug mode) 
%
%
%   Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA 2009.

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.4                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('----------------------- Pessoa: Target Set Initiated ------------------- ');

nbatch=N_sims;

load(strcat(system_name, '_symb'),'params_symb');

% transform target set to labels

% Target Set has to be Under-approximated
tmin=zeros(params_symb.n,1);
tmax=zeros(params_symb.n,1);
for k=1:params_symb.n
    if (mod(target_set(k,1),params_symb.eta)~=0)
        tmin(k)=ceil(target_set(k,1)/params_symb.eta);
    else
        tmin(k)=target_set(k,1)/params_symb.eta;
    end

    if (mod(target_set(k,2),params_symb.eta)~=0)
        tmax(k)=floor(target_set(k,2)/params_symb.eta);
    else
        tmax(k)=target_set(k,2)/params_symb.eta;
    end
end
maxc=tmax-params_symb.min(params_symb.xoind);
minc=tmin-params_symb.min(params_symb.xoind);

%% Check here for max and min >0 and <params_symb.num

if( sum(maxc<0) || sum(minc<0) || sum(maxc>params_symb.num(params_symb.xoind)) || sum(minc>params_symb.num(params_symb.xoind)))
   error('The target set must be a subset of the state set.'); 
end

%%
totloops=ceil(prod((maxc-minc)+ones(size(maxc)))/nbatch);

% Here we call the pss_target_set to find the target set
if(type)
    pessoa_charf2bdd(filename_set,minc,maxc,verbose);
else
    pessoa_set2bdd(filename_set,minc,maxc,verbose);
end

filename=strcat(filename_set,'_symb');
save(filename, 'params_symb');
