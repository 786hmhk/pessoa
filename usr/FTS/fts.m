%% FTS Example
clear all;
close all;

%% System definition

FTS_cont.custom=1;

%% Abstraction parameters 

FTS_symb.xset = [0 4];   % 5 states
FTS_symb.uset = [0 1];   % One fault
                                                

FTS_symb.tau = 1; % Time quantization
FTS_symb.eta = 1;   % State quantization
FTS_symb.mu = 1;    % Input quantization


%% Abstract


FTS_name='FTS';

Relation=1;

BatchSize=1;

Verbose=1;

% The command pss_abstract computes the symbolic model.

pss_abstract(FTS_cont,FTS_symb,FTS_name,Relation,BatchSize,Verbose);

%% Controller design

TargetSet_name{1}='p'; % SA
%
%
TargetSet=[4 4];

TargetSetType=0;
%
% % Build the target set
%
pss_build_set(FTS_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)


FTSController_name='FTSController';

ControlObj = 2;

pss_design(FTS_name,TargetSet_name,FTSController_name,ControlObj,BatchSize,Verbose)
%% Print traces

for(i=0:FTS_symb.xset(2))
   display('traces from initial state');
   display(i);
   iterate_fti(FTSController_name,i,[],1,1)
end


%% Abstract


FTS_name='FTS_nd';

Relation=2;

BatchSize=1;

Verbose=1;

% The command pss_abstract computes the symbolic model.

pss_abstract(FTS_cont,FTS_symb,FTS_name,Relation,BatchSize,Verbose);


%% Controller design

TargetSet_name{1}='p'; % SA
%
%
TargetSet=[4 4];

TargetSetType=0;
%
% % Build the target set
%
pss_build_set(FTS_name,TargetSet,TargetSet_name{1},BatchSize,TargetSetType,Verbose)


FTSController_name='FTSController_nd';

ControlObj = 2;

pss_design(FTS_name,TargetSet_name,FTSController_name,ControlObj,BatchSize,Verbose)

%% Print traces

for(i=0:FTS_symb.xset(2))
   display('traces from initial state');
   display(i);
   iterate_fti(FTSController_name,i,[],1,1)
end
