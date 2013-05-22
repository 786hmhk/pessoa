clear all;
close all;

%% SYSTEM DEFINITION

system1_cont.custom=1;

%% ABSTRACTION PARAMETERS

system1_symb.tau=0.5;
system1_symb.eta=0.2;
system1_symb.mu=0.5;
system1_symb.xset=[1 5; 1 5; -pi-system1_symb.eta pi+system1_symb.eta];
system1_symb.uset=[0 0.5; -0.5 0.5];


%% ABSTRACT %%

% Verbose Level 0 nothing info, 1 messages printed on screen and profiler, 2 as 1 but also plots resulting FSMs (just works for 2D).
verbose=1;
% Size of simulation batches
N_sims=3000;
% Label 1 for deterministic system, 2 for non-deterministic
deter=2;
% system name
system1_name='unicycle';

% Construct abstraction
pss_abstract(system1_cont,system1_symb,system1_name,deter,N_sims,verbose)

%% DESIGN %%

% Verbose Level 0 nothing info, 1 messages printed on screen and profiler, 2 as 1 but also plots resulting FSMs (just works for 2D).
verbose=1;
% Size of simulation batches
N_sims=1000;
% system name
system1_name='unicycle';

% Target set now defined through a characteristic function (pss_target_set)
% Set to entire state space

% Desired target_set = {x\in [4.6,5], y\in[1,1.6]} 
target_set = [4.6 5; 1 1.6; -pi-system1_symb.eta pi+system1_symb.eta];

% Build the specification system (a set in this case).
% Construct_type: 0 - using hyperboxes, 1 - using pss_target_set.m
construct_type=0;
% Set filename
set_nameT='Target';

% %Build Target Set
pss_build_set(system1_name,target_set,set_nameT,N_sims,construct_type,verbose)

% Constraint Set, i.e. definition of obstacles, via a characteristic
% function
construct_type=1;
set_nameC='Constraint';
pss_build_set(system1_name,system1_symb.xset,set_nameC,N_sims,construct_type,verbose)


set_name{1}=set_nameT;
set_name{2}=set_nameC;
% Design the controller
contr_name='unicyle_control';

 type=4; % Reach and Stay Target Set while Staying in Safe Set
 pss_design(system1_name,set_name,contr_name,type,N_sims,verbose)


%% TEST CONTROL %%

xini=[1.4;1;0];
uini=[0;0];

contr_name='unicyle_control';
Robot_smlnk;
