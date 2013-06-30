

clear all;


%shortestPath();







%% DC Motor Example
%  Paulo Tabuada 
%  CyPhyLab, October 2009

clear all;
close all;

%% System definition

% Physical parameters

J=0.00025;
B=0.0001;
k=0.05;
L=0.0015;
R=0.5;

% Continuous dynamics (linear)

DCMotor_cont.A = [-B/J k/J; -k/L -R/L];     % A matrix
DCMotor_cont.B=[0; 1/L];                    % B matrix

%% Abstraction parameters

DCMotor_symb.xset=[-1 30; -10 10];         % Range of values for the states
                                           % in the symbolic model.
                                           
DCMotor_symb.uset=[-10 10];                % Range of values for the input 
                                           % in the symbolic model.
                                           
DCMotor_symb.tau=0.05;                     % Time quantization
DCMotor_symb.eta=0.5;                      % State quantization
DCMotor_symb.mu=0.01;                      % Input quantization


%% Abstract

% The variable DCMotor_name records the name of the file where the
% symbolic model will be stored.

DCMotor_name='DCMotor';

% The Relation flag controls what kind of relation is to be established
% between the original model and the symbolic model:
%       Relation=1 creates a symbolic model related to the original model
% by an approximate bisimulation relation.
%       Relation=2 creates a symbolic model related to the original model
% by an approximate alternating simulation relation.

% For the chosen values of tau, eta, and mu we have an approximate
% bisimulation. Hence we set Relation to 1.

Relation=1;

% The BatchSize flag controls the size of the computation batches.
% Transitions are generated and added to the symbolic model in batches of
% size BatchSize. Different values for BatchSize affect the computation
% time of the symbolic model. Typical values range between 100 and 1000.

BatchSize=1000;

% The Verbose flag lets you control how much information is displayed while
% the symbolic model is being computed:
%       Verbose=0 provides no information.
%       Verbose=1 provides messages printed on the screen.
%       Verbose=2 provides a graphical representation for the symbolic
% model in addition to the messages printed when the flag is set to 1. 
% This flag is only available for 2-dimensional systems.
%       Verbose=3 debug mode

Verbose=2;

% The command pss_abstract computes the symbolic model.

%pss_abstract(DCMotor_cont,DCMotor_symb,DCMotor_name,Relation,BatchSize,Verbose)

%% Synthesize controller

% Verbose flag. Same as before.
Verbose=3;

% BatchSize flag. Same as before.
BatchSize=1000;

SystemCost_name = 'DCMotorCosts';

% Create the systems cost ADD
pss_build_cost(DCMotor_name, SystemCost_name, BatchSize, Verbose);






