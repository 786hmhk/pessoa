function pss_abstract(name_cont,name_symb,filename_symb_sys,deter,N_sims,verbose)

%
% pss_abstract_add(system_cont,system_symb,filename,relation,batch_size,verbose)
%
% Constructs a finite abstraction of a continuous-time linear control system
% and stores it as a reduced ordered binary decision diagram.
%   
% INPUTS: system_cont  - structure describing the control system (see online
%                        documentation);
%         system_symb  - structure describing the abstraction parameters (see 
%                        online documentation);
%         filename     - string containing the name of the file where the 
%                        abstraction will be saved (2 files will be created,
%                        filename.bdd and filename_symb.mat);
%         relation     - flag defining the type of abstraction: 
%                        relation=1 requests the approximate alternating 
%                                   bisimulation algorithm;
%                        relation=2 requests the approximate alternating 
%                                   simulation algorithm;
%         batch_size   - size of computation batches to be used (see online
%                        documentation);
%         verbose      - flag defining the verbosity level: 
%                        verbose=0 minimal;
%                        verbose=1 extra information;
%                        verbose=2 several plots in addition to the
%                                  information given by verbose=1; (only 
%                                  works for 2D systems);
%                        verbose=3 debug mode.
%
% OUTPUT: none
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab May 2010.

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.4                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('----------------------- Pessoa: Abstraction Initiated ------------------');

%name_symb.xset,name_symb.uset,name_symb.tau,name_symb.eta,name_symb.mu,name_symb.deter;
%name_cont.A,name_cont.B,name_cont.C,name_cont.xset,name_cont.uset,name_cont.dset;

if ~isfield(name_symb,'epsilon')
    name_symb.epsilon=name_symb.eta/2;
end

if isfield(name_cont,'num_sw')
    name_symb.uset=[name_symb.uset;name_symb.mu*[1 name_cont.num_sw]];
    params_symb=pss_build_params_symb(name_symb.xset,name_symb.uset,name_symb.tau,name_symb.eta,name_symb.mu,name_symb.epsilon,deter,N_sims);
    params_symb.num_sw=name_cont.num_sw;
else
    params_symb=pss_build_params_symb(name_symb.xset,name_symb.uset,name_symb.tau,name_symb.eta,name_symb.mu,name_symb.epsilon,deter,N_sims);
end;

% If the system is linear and with disturbances
if isfield(name_cont,'A') && isfield(name_cont,'B') && isfield(name_symb,'dset') && isfield(name_cont,'C')
    params_cont=pss_build_params_cont(params_symb,name_cont.A,name_cont.B,name_cont.C,name_symb.dset);
elseif isfield(name_cont,'A') && isfield(name_cont,'B') %If it is only linear
    params_cont=pss_build_params_cont(params_symb,name_cont.A,name_cont.B);
end

if isfield(name_cont,'custom')
    params_cont.custom=name_cont.custom;
else
    params_cont.custom=0;
end

nbatch=N_sims;
totloops=ceil(prod(params_symb.num+ones(size(params_symb.num)))/nbatch);

filename=strcat(filename_symb_sys,'_symb');
save(filename, 'params_symb');

if isfield(name_cont,'A') && isfield(name_cont,'B') && ~isfield(name_cont,'num_sw') %If it is linear (non switched system)
	xoob=zeros(params_symb.n,1); xoob(1)=params_symb.num(1)+1;
	Ad = params_cont.Ad;
	Bd = params_cont.Bd;
	minp = params_cont.minp;
	vmax = params_cont.vmax;
	%vmin = params_cont.vmin;
	xoffset = params_symb.min(params_symb.xoind);
	uoffset = params_symb.min(params_symb.uind);
	xset = name_symb.xset;
	uset = name_symb.uset;
	tau = name_symb.tau;
	mu	= name_symb.mu;
	eta = name_symb.eta;

	pessoa_abstract_add(filename_symb_sys,verbose,Ad,Bd,minp,vmax,xoffset,uoffset,xoob,xset,uset,tau,mu,eta); 
else
	pessoa_abstract(filename_symb_sys,verbose);
end


