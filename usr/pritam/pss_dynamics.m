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

%% INITIALIZATIONS

% lmaxtemp=zeros(params_symb.n,1);
% lmintemp=lmaxtemp;
% 
% %% SIMULATION
% 
% 
% lmaxtemp_row=params_cont.Ad*xl+params_cont.Bd*ul+params_cont.minp+params_cont.vmax;
% lmintemp_row=params_cont.Ad*xl+params_cont.Bd*ul+params_cont.minp+params_cont.vmin;
% 
% % This is the check to see if the intersection with neighbouring cells
% % is empty or not.
% for k=1:params_symb.n
%     if(mod(lmaxtemp_row(k),1)~=1/2)
%         lmaxtemp(k)=round(lmaxtemp_row(k));
%     else
%         lmaxtemp(k)=floor(lmaxtemp_row(k));
%     end;
%     if(mod(lmintemp_row(k),1)~=1/2)
%         lmintemp(k)=round(lmintemp_row(k));
%     else
%         lmintemp(k)=ceil(lmintemp_row(k));
%     end;
% end;
% % Shift to get labels between 0 and maximum label
% lmax=lmaxtemp-params_symb.min(params_symb.xoind);
% lmin=lmintemp-params_symb.min(params_symb.xoind);
% %%%
% 
% %%% Add all the transitions generated:
% x=lmin;
% ntrans=1;
% jl=1;
% while(jl<=params_symb.n && x(params_symb.n)<=lmax(params_symb.n))
% jl=1;
%     xendtemp(:,ntrans) = x;
%     while(jl<params_symb.n && x(jl)>=lmax(jl))
%         x(jl)=lmin(jl);
%         jl=jl+1;
%     end
%     x(jl)=x(jl)+1;
%     ntrans=ntrans+1;
% end
% xend=xendtemp(:,1:ntrans-1);

%% All by hand

eta=params_symb.eta;
etah=eta/2;

x = (xl+params_symb.min(params_symb.xoind))*params_symb.eta;
u = (ul+params_symb.min(params_symb.uind))*params_symb.mu;

% disturbance limits
udmax=0.1;
udmin=-0.1;

x1max= x(1)+etah +(x(3)-x(2)+eta)*params_symb.tau + (udmax-u)*params_symb.tau^2/2;
x2max= x(2)+etah +u*params_symb.tau;
x3max= x(3)+etah +udmax*params_symb.tau;

x1min= x(1)-etah +(x(3)-x(2)-eta)*params_symb.tau + (udmin-u)*params_symb.tau^2/2;
x2min= x(2)-etah +u*params_symb.tau;
x3min= x(3)-etah +udmin*params_symb.tau;

% x1max= x(1)+etah +(x(2)+eta)*params_symb.tau;
% x2max= x(2)+etah +(udmax-u)*params_symb.tau;
% x3max= x(3)+etah +u*params_symb.tau;
% 
% x1min= x(1)-etah +(x(2)-eta)*params_symb.tau;
% x2min= x(2)-etah +(udmin-u)*params_symb.tau;
% x3min= x(3)-etah +u*params_symb.tau;

xmax=[x1max;x2max;x3max];
xmin=[x1min;x2min;x3min];

xmax_label=zeros(params_symb.n-params_symb.m,1);
xmin_label=xmax_label;

%compute the labels for max and min
for k=1:params_symb.n
    if (mod(xmax(k),eta)~=etah)
        xmax_label(k)=round(xmax(k)/eta)-params_symb.min(k);
    else
        xmax_label(k)=floor(xmax(k)/eta)-params_symb.min(k);
    end

    if (mod(xmin(k),eta)~=etah)
        xmin_label(k)=round(xmin(k)/eta)-params_symb.min(k);
    else
        xmin_label(k)=ceil(xmin(k)/eta)-params_symb.min(k);
    end
end

% attach transitions
xend=[];
x=xmin_label;
ntrans=1;
jl=1;
while(jl<=params_symb.n && x(params_symb.n)<=xmax_label(params_symb.n))
jl=1;
    xend=[xend x];
    while(jl<params_symb.n && x(jl)>=xmax_label(jl))
        x(jl)=xmin_label(jl);
        jl=jl+1;
    end
    x(jl)=x(jl)+1;
    ntrans=ntrans+1;
end

% saturate at max and min velocities the target vehicle and the
% intervehicle distance measurements
 
xend(1,xend(1,:)>params_symb.num(1))=params_symb.num(1);
%xend(1,xend(1,:)<0)=0;
%xend(2,xend(2,:)>params_symb.num(2))=params_symb.num(2);
%xend(2,xend(2,:)<0)=0;
xend(3,xend(3,:)>params_symb.num(3))=params_symb.num(3);
xend(3,xend(3,:)<0)=0;

