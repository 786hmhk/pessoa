%% DC-DC converter Example
%  Manuel Mazo Jr.
%  CyPhyLab, May 2010

clear all;
close all;

%% System definition

% Physical parameters

vs=1;
L=3;
C=70;
RL=0.05;
R0=1;
RC=0.005;

% Continuous dynamics (linear)

DCDC_cont.A = {[-RL/L 0; 0 -(1/C)*(1/(RC+R0))], [-(1/L)*(RL+(RC*R0)/(RC+R0)) -(1/5)*(1/L)*R0/(RC+R0); 5*(1/C)*R0/(RC+R0) -(1/C)/(RC+R0)]};     % A matrices
DCDC_cont.B={[1/L;0], [1/L;0]};                    % B matrices
DCDC_cont.num_sw=2;

%% Abstraction parameters

DCDC_symb.xset=[1.2 1.6; 5.6 5.8];         % Range of values for the states
                                           % in the symbolic model.
                                           
DCDC_symb.uset=[vs vs];                % Range of values for the input 
                                           % in the symbolic model.

DCDC_symb.tau=0.2;                     % Time quantization
DCDC_symb.eta=0.0071;                      % State quantization
DCDC_symb.mu=1;                      % Input quantization
DCDC_symb.epsilon=3;

% DCDC_symb.tau=0.5;
% DCDC_symb.eta=7e-4;
% DCDC_symb.mu=1;
% DCDC_symb.epsilon=0.5;


%% Abstract

% The variable DCDC_name records the name of the file where the
% symbolic model will be stored.

DCDC_name='DCDC';

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

pss_abstract(DCDC_cont,DCDC_symb,DCDC_name,Relation,BatchSize,Verbose)

%% Synthesize controller

% Verbose flag. Same as before.
Verbose=1;

% BatchSize flag. Same as before.
BatchSize=1000;

% The parameter TargetSet_name records the name of the file where the
% target set will be stored.

TargetSet_name{1}='DCDCTargetSet';

% The target set can be defined by a rectangle (faster) or by defining the
% characteristic function of the desired target set (slower). The parameter
% TargetSet is used to specify the target set rectangle. To define the
% target set through its characteristic function, please edit the file 
% pss_target_set. When using the characteristic function method, it is
% assumed that the parameter TargetSet is an over-approximation of the 
% the target set.

TargetSet=[1.2 1.6; 5.6 5.8];

% The flag TargetSetType determines if a rectangle or a characteristic
% function is to be used:
%       TargetSetType=0: rectangle.
%       TargetSetType=1: characteristic function.

TargetSetType=0;

% Build the target set

pss_build_set(DCDC_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)

% The parameter DCDCController_name records the name of the file where
% the symbolic controller will be stored.

DCDCController_name='DCDCController';

% The flag ControlObjective specifies the control objective as a function
% of two auxiliary sets: the target set W and the constraint set Z.
%       ControlObjective=1: remain within the target set (always W).
%       ControlObjective=2: reach the target set (eventually W).
%       ControlObjective=3: reach and remain within the target set
% (eventually always W).
%       ControlObjective=4: reach and remain within the target set while 
% always remaining within the constraint set (eventually always W and 
% always Z).

ControlObjective=1;

% Synthesize the controller

pss_design(DCDC_name,TargetSet_name,DCDCController_name,ControlObjective,BatchSize,Verbose)

%% Simulate closed-loop

A1=DCDC_cont.A{1};
B1=DCDC_cont.B{1};
A2=DCDC_cont.A{2};
B2=DCDC_cont.B{2};
A={A1,A2};
B={B1,B2};
C=eye(size(A1));
D=0*B1;
xini=[1.35;5.65];
uini=[1;1];
DCDCSmlnk;
