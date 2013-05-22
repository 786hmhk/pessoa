clear all;
close all;

%% SYSTEM DEFINITION

system1_cont.custom=1;

%% ABSTRACTION PARAMETERS

system1_symb.tau=0.5;
system1_symb.eta=0.2;
system1_symb.mu=0.2;
system1_symb.xset=[1 5; 1 5; -pi-system1_symb.eta pi+system1_symb.eta; 0 0.4; -0.4 0.4];
system1_symb.uset=[0 0.4; -0.4 0.4];


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
