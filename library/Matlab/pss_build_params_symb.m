function params_symb=pss_build_params_symb(xset,uset,tau,eta,mu,epsilon,deter,N_sims)

%
%params_symb=pss_build_params_symb(xset,uset,tau,eta,mu,deter,N_sims)
%
%Fills the structure "params_symb" containing information about the symbolic abstraction
%to be constructed, and used in the construction.
%   
%INPUTS: xset	- matrix defining the state set, each row contains the minimum and maximum value for that coordinate.
%        uset	- matrix defining the input set, each row contains the minimum and maximum value for that coordinate.
%        tau	- sampling time used in the abstraction.
%        eta	- space sampling unit used in the abstraction.
%		 mu		- input sampling unit used in the abstraction.
%		 deter	- flag to select deterministic (1) or non-deterministic (2) construction.
%		 N_sims	- size of batches to be used in the abstraction construction.
%
%OUTPUTS: params_symb	- resulting structure.
%
%For more details see the Online Documentation.
%
%Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA May 2010.


% params_symb.n
% params_symb.m
% params_symb.tau
% params_symb.eta
% params_symb.mu
% params_symb.xoind
% params_symb.uind
% params_symb.xeind
% params_symb.min
% params_symb.max
% params_symb.num
% params_symb.nume
% params_symb.nxbits
% params_symb.nubits
% params_symb.nbitsx
% params_symb.nbitsloop
% params_symb.totbits
% params_symb.nbits
% params_symb.deter


% Semantics

params_symb.n=size(xset,1);		% dim(x)
params_symb.m=size(uset,1);		% dim(u)
params_symb.tau=tau;			% tau
params_symb.eta=eta;			% nu
params_symb.mu=mu;				% mu
params_symb.epsilon=epsilon;    % epsilon

% If bisimulation indicated, check if the parameters satisfy it

                         
% Codification of states/inputs, numbits, etc.

params_symb.xoind=[1:params_symb.n];
params_symb.uind=[params_symb.n+1:params_symb.n+params_symb.m];
params_symb.xeind=[params_symb.n+params_symb.m+1:2*params_symb.n+params_symb.m];

% Input Set is Under-approximated
umin=zeros(params_symb.m,1);
umax=zeros(params_symb.m,1);
for k=1:params_symb.m
    if (mod(uset(k,1),params_symb.mu)~=0)
        umin(k)=ceil(uset(k,1)/params_symb.mu);
    else
        umin(k)=uset(k,1)/params_symb.mu;
    end

    if (mod(uset(k,2),params_symb.mu)~=0)
        umax(k)=floor(uset(k,2)/params_symb.mu);
    else
        umax(k)=uset(k,2)/params_symb.mu;
    end
end
% State Set is Over-approximated
xmin=zeros(params_symb.n,1);
xmax=zeros(params_symb.n,1);
for k=1:params_symb.n
    if (mod(xset(k,1),params_symb.eta)~=0)
        xmin(k)=floor(xset(k,1)/params_symb.eta);
    else
        xmin(k)=xset(k,1)/params_symb.eta;
    end

    if (mod(xset(k,2),params_symb.eta)~=0)
        xmax(k)=ceil(xset(k,2)/params_symb.eta);
    else
        xmax(k)=xset(k,2)/params_symb.eta;
    end
end

params_symb.min=round([xmin; umin]);
params_symb.max=round([xmax; umax]);

params_symb.num=params_symb.max-params_symb.min;
params_symb.nume=[params_symb.num; params_symb.num(params_symb.xoind)+[1;zeros(params_symb.n-1,1)]]; 
% Line above modified so that we'll plot out-of-bounds transitions, see below.

   params_symb.nxbits(1)=length(dec2bin(1+params_symb.num(1))); % For the out-of-bounds state, namely oobs=[max(1)+1; 0;...;0]
for i=2:params_symb.n
   params_symb.nxbits(i)=length(dec2bin(params_symb.num(i)));
end
for i=params_symb.n+1:params_symb.n+params_symb.m
   params_symb.nubits(i-params_symb.n)=length(dec2bin(params_symb.num(i)));
end

params_symb.nbitsx=sum([params_symb.nxbits]);
params_symb.nbitsloop=sum([params_symb.nxbits params_symb.nubits]);

params_symb.nbits=[params_symb.nxbits params_symb.nubits params_symb.nxbits];
params_symb.totbits=sum([params_symb.nxbits params_symb.nubits params_symb.nxbits]);

params_symb.deter=deter;
