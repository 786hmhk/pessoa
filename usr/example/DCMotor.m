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

Verbose=1;

% The command pss_abstract computes the symbolic model.

pss_abstract_add(DCMotor_cont,DCMotor_symb,DCMotor_name,Relation,BatchSize,Verbose)

%% Synthesize controller

% Verbose flag. Same as before.
Verbose=1;

% BatchSize flag. Same as before.
BatchSize=1000;

% The parameter TargetSet_name records the name of the file where the
% target set will be stored.

TargetSet_name{1}='DCMotorTargetSet';

% The target set can be defined by a rectangle (faster) or by defining the
% characteristic function of the desired target set (slower). The parameter
% TargetSet is used to specify the target set rectangle. To define the
% target set through its characteristic function, please edit the file 
% pss_target_set. When using the characteristic function method, it is
% assumed that the parameter TargetSet is an over-approximation of the 
% the target set.

TargetSet=[19.5 20.5; -10 10];

% The flag TargetSetType determines if a rectangle or a characteristic
% function is to be used:
%       TargetSetType=0: rectangle.
%       TargetSetType=1: characteristic function.

TargetSetType=0;

% Build the target set

pss_build_set(DCMotor_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)

% The parameter DCMotorController_name records the name of the file where
% the symbolic controller will be stored.

DCMotorController_name='DCMotorController';

% The flag ControlObjective specifies the control objective as a function
% of two auxiliary sets: the target set W and the constraint set Z.
%       ControlObjective=1: remain within the target set (always W).
%       ControlObjective=2: reach the target set (eventually W).
%       ControlObjective=3: reach and remain within the target set
% (eventually always W).
%       ControlObjective=4: reach and remain within the target set while 
% always remaining within the constraint set (eventually always W and 
% always Z).

ControlObjective=3;

% Synthesize the controller

pss_design(DCMotor_name,TargetSet_name,DCMotorController_name,ControlObjective,BatchSize,Verbose)

%% Simulate closed-loop

A=DCMotor_cont.A;
B=DCMotor_cont.B;
C=eye(size(A));
D=0*B;
xini=[0;0];
uini=0;
DCMotorSmlnk;

%% New model for the DC motor circuit

% We no longer treat the voltage as the input variable but instead assume
% the DC motor circuit to be connected to a 10V source through an H-bridge.
% We control directly the opening and closing of the bridge switches.

% We redefine the quantization parameters. The time quantization will have
% to be smaller to capture the very frequent operation of the switches.
% This requires the space quantization to also be reduced. The input
% quantization is increased to 10 so that only three voltage values can be
% applied: -10, 0, and 10. These values correspond to the several switch
% configurations in the H-bridge.

DCMotor_symb.tau=0.0001;                     
DCMotor_symb.eta=0.05;                       
DCMotor_symb.mu=10;                          

% For the new quantization parameters we can only obtain an approximate 
% alternating simulation relation. Hence, we change the flag Relation.

Relation=2;

% We save the new abstraction in a different file.

DCMotor_name='DCMotorSW';

% And compute the abstraction.

pss_abstract(DCMotor_cont,DCMotor_symb,DCMotor_name,Relation,BatchSize,Verbose)

%% Synthesize a new controller

% We construct a new target set since the new abstraction has a different
% state quantization (eta). The new target set will be saved in a different
% file.

TargetSet_name{1}='DCMotorSWTargetSet';

% Build the target set

pss_build_set(DCMotor_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)

% We save the new controller on a different file.

DCMotorController_name='DCMotorSWController';

% Synthesize the controller

pss_design(DCMotor_name,TargetSet_name,DCMotorController_name,ControlObjective,BatchSize,Verbose)

%% Simulate new closed-loop

DCMotorSmlnk;

%% Synthesize a controller with state constraints

% We construct a new target set where the current is restricted to range
% between -2 and 2 Amperes.

TargetSet=[19.5 20.5; -0.7 0.7];
TargetSet_name{1}='DCMotorSWTargetSet2';

% Build the target set

pss_build_set(DCMotor_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)

% We construct the constraint set that will impose state constraints 
% holding for all the time.

TargetSet_name{2}='DCMotorSWConstraintSet';

TargetSet=[-1 30;-3 3]; 

pss_build_set(DCMotor_name,TargetSet,TargetSet_name{2},BatchSize,TargetSetType,Verbose)

% The new requirement will be of type 4: eventually always W and always Z.

ControlObjective=4;

% We save the new controller on a different file.

DCMotorController_name='DCMotorSWController2';

% Synthesize the controller

pss_design(DCMotor_name,TargetSet_name,DCMotorController_name,ControlObjective,BatchSize,Verbose)

%% Simulate new closed-loop

DCMotorSmlnk;

