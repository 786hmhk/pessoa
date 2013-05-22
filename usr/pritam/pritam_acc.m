%% ACC Example
clear all;
close all;

%% System definition

% Physical parameters

% Continuous dynamics (linear)
% assume host vehicle is following a target vehicle

% state vector sv = [xr vr vh]'
% xr = relative distance between host and target vehicle  (in meter)
% vr = relative velocity between host and target vehicle  (in m/s)
% vh = velocity of the host                       (in m/s)
% vt = velocity of the target                          (in m/s)
% acch = accl of host vehicle                               (in m/s.s)
% acct = accl of target vehicle (unknown)                   (in m/s.s)

ACC_cont.A = [0 -1 1;0 0 0; 0 0 0];    % A matrix
ACC_cont.B=[0;1;0]; % B matrix
ACC_cont.C = [0;0;1]; % C matrix (for disturbance input)

ACC_cont.custom=1;

%% Abstraction parameters 

ACC_symb.xset = [0 70; 25 35; 28 30];   % Range of values for the states
                                                % in the symbolic model.
ACC_symb.uset = [-0.4,0.4];                         % Range of values for the input (host acceleration)
                                                % in the symbolic model.
                                                
% The disturbance range is set in pss_dynamics                                                
ACC_symb.dset = [-0.1,0.1];                     % Range of values for the disturbance (target acceleration)
                                                % in the symbolic model.
ACC_symb.tau = 0.00001; % Time quantization
ACC_symb.eta = 0.5;   % State quantization
ACC_symb.mu = 0.1;    % Input quantization


%% Abstract


ACC_name='ACC';

Relation=2;

BatchSize=10000;

Verbose=1;

% The command pss_abstract computes the symbolic model.

pss_abstract(ACC_cont,ACC_symb,ACC_name,Relation,BatchSize,Verbose);

%% Controller design
TargetSet_name{1}='p'; % SA
%
%
TargetSet=[30 70; 27 35; 28 30];

TargetSetType=0;
%
% % Build the target set
%
pss_build_set(ACC_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)

ACCController_name='ACCController';

ControlObj = 1;

pss_design(ACC_name,TargetSet_name,ACCController_name,ControlObj,BatchSize,Verbose)

%% SIMULATE!! :)

A=ACC_cont.A;
B=ACC_cont.B;
Bd=ACC_cont.C;
C=eye(size(A));
D=0*B;
xini=[40;29;30];
uini=0;
ACC;
