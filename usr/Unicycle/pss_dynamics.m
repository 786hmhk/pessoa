function xend=pss_dynamics(params_symb, params_cont, xl, ul)

%
%   xend = pss_dynamics(params_symb, params_cont, xl, ul)
%
%   This function computes the transitions in labels (as opposed to coordinates) directly. 
%	In general after the simulation step one would need to call "pss_build_fsm_indices.m"
%   at the end to return labels.
%   
%       INPUTS: xl - state labels.
%           	ul - input labels.
%
%       OUTPUT: xend	- ending state(s) in labels.
%
%   Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA 2009.

%% Transform labels into states (I'll probably make a function for this like the one at the end)

u=params_symb.mu*(ul+params_symb.min(params_symb.uind)); %u(1)=v, u(2)=w
x=params_symb.eta*(xl(1:params_symb.n-params_symb.m)+params_symb.min(1:params_symb.n-params_symb.m)); %x(1)=x, x(2)=y, x(3)=\theta

%% Here computations on states

% Some previous computations and other arrangements neccessary
minpoint_label=zeros(params_symb.n-params_symb.m,1);
maxpoint_label=minpoint_label;
tau=params_symb.tau;
eta=params_symb.eta;
etah=params_symb.eta/2;

% computation of the reachable set

wt=u(2)*tau;
wth=wt/2;
tpwth=x(3)+wth;


if u(2)~=0
    kvw=2*(u(1)/u(2))*sin(wth);    
else    
    kvw=u(1)*tau;
end

x1plus1=x(1)+kvw*cos(x(3)+wth+etah);
x1plus2=x(1)+kvw*cos(x(3)+wth-etah);

x2plus1=x(2)+kvw*sin(x(3)+wth+etah);
x2plus2=x(2)+kvw*sin(x(3)+wth-etah);

%possible extrema inside the interval

if abs(x(3)+wth)<etah
    x1plus3=x(1)+kvw;    
elseif abs(x(3)+wth-pi)<etah || abs(x(3)+wth+pi)<etah
    x1plus3=x(1)-kvw;    
else
    x1plus3=x1plus2;
end

if abs(x(3)+wth-pi/2)<etah
    x2plus3=x(2)+kvw;
elseif abs(x(3)+wth+pi/2)<etah
    x2plus3=x(2)-kvw;
else
    x2plus3=x2plus2;
end

x1plusMax=max([x1plus1,x1plus2,x1plus3])+etah;
x1plusMin=min([x1plus1,x1plus2,x1plus3])-etah;
x2plusMax=max([x2plus1,x2plus2,x2plus3])+etah;
x2plusMin=min([x2plus1,x2plus2,x2plus3])-etah;
x3plusMax=x(3)+wt+etah;
x3plusMin=x(3)+wt-etah;

maxpoint=[x1plusMax; x2plusMax; x3plusMax];
minpoint=[x1plusMin; x2plusMin; x3plusMin];

% Label computations

for k=1:params_symb.n-params_symb.m
    if (mod(maxpoint(k),eta)~=etah)
        maxpoint_label(k)=round(maxpoint(k)/eta)-params_symb.min(k);
    else
        maxpoint_label(k)=floor(maxpoint(k)/eta)-params_symb.min(k);
    end

    if (mod(minpoint(k),eta)~=etah)
        minpoint_label(k)=round(minpoint(k)/eta)-params_symb.min(k);
    else
        minpoint_label(k)=ceil(minpoint(k)/eta)-params_symb.min(k);
    end
end

% Add all the points in the Post
% Some check needed to wrap around in the 3rd coordinate

xend=[];
x=minpoint_label;
ntrans=1;
jl=1;
while(jl<=params_symb.n-params_symb.m && x(params_symb.n-params_symb.m)<=maxpoint_label(params_symb.n-params_symb.m))
jl=1;
    xend=[xend x];
    while(jl<params_symb.n-params_symb.m && x(jl)>=maxpoint_label(jl))
        x(jl)=minpoint_label(jl);
        jl=jl+1;
    end
    x(jl)=x(jl)+1;
    ntrans=ntrans+1;
end


% Attach the output of the second state (input at previous state)

xend=[xend; kron(ul,ones(1,size(xend,2)))];

%% Transformation back to labels
% Transform the vector xend into labels (I'll make it soon so that
% pss_build_fsm_indices operates directly on vectors so that this loop is
% inside it

% Not needed we computed labels already!!
% for i=1:size(xend,2)
%      xend(:,i)=pss_build_fsm_indices(params_symb,xend(:,i),'x');
% end
