function [trans xnext unext]=pss_sim_transit(xini,u,N,params_symb,params_cont)

%
% bddstring=pss_sim_transit(xini,u,batch_size,params_symb,params_cont)
%
% Simulator of the system dynamics. From the discretized dynamics of the
% linear system with time step "tau", this function computes the 
% end-point(s) ("xend", internal variable) given initial states "xini" and
% inputs "u". It computes these transitions for "batch_size" combinations
% "xini,u" starting from the ones provided. This function also creates all
% the transitions to be added to the symbolic model abstraction, i.e. 
% (xini)x(u)x(xend_i). The provided state and input, and the generated
% output (transitions), are represented as labels and not as real coordinates.
%   
% INPUTS: xini       - initial state represented in labels;
%         u          - input represented in labels;
%         batch_size - size of computation batches to be used (see online
%						documentation);
%         params_symb - structure containing the abstraction parameters
%                       (see online documentation);
%         params_cont - structure containing the description of the 
%                       continuous-time model (see online documentation).
%
% OUTPUT: trans - matrix containing the column vectors [xo u xe]' 
%                 representing the transitions to be added to the
%                 abstraction under construction;
%         xnext - initial state to be used in the next batch computation
%                 (in labels);
%         unext - input to be used in the next batch computation (in labels).
%
% CALLS:  pss_dynamics.m
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab May 2010.

%% INITIALIZATIONS
npm=params_symb.n+params_symb.m;
xuo=[xini;u];
xoob=zeros(params_symb.n,1); xoob(1)=params_symb.num(1)+1;
trans_tmp=[];
xendtemp=[];

if isfield(params_cont,'Ad') && isfield(params_cont,'Bd') && ~params_cont.custom %LINEAR SYSTEMS NOT-CUSTOM
    newrow=1;
    %% DETERMINISTIC LINEAR SYSTEMS
    if params_symb.deter==1
        trans=[];
    else %% NON-DETERMINISTIC LINEAR SYSTEMS
%        trans_tmp=zeros(npm+params_symb.n,N*prod(ceil((params_cont.vmax-params_cont.vmin)./params_symb.eta)+1));
%        xendtemp=zeros(params_symb.n,prod(ceil((params_cont.vmax-params_cont.vmin)./params_symb.eta)+1));
        lmaxtemp=zeros(params_symb.n,1);
        lmintemp=lmaxtemp;
        counter=1;
    end;
else
    trans=[];
end

%% SIMULATIONS
for i=1:N
% The simulation of the system is performed in the next lines, the result is also
% transformed into labels in N+.

xl=xuo(params_symb.xoind);
ul=xuo(params_symb.uind);

    if isfield(params_cont,'Ad') && isfield(params_cont,'Bd') && ~params_cont.custom %LINEAR SYSTEMS NOT-CUSTOM
    
        %% DETERMINISTIC LINEAR SYSTEMS
        if params_symb.deter==1

            if isfield(params_cont,'num_sw')
                if newrow
                    sim_row=(params_cont.Ad{ul(end)+1}*xl+params_cont.Bd{ul(end)+1}*ul(1:end-1)+params_cont.minp{ul(end)+1});
                    xend=round(sim_row)-params_symb.min(params_symb.xoind);
                    newrow=0;
                else
                    sim_row=sim_row+params_cont.Ad{ul(end)+1}(:,1);
                    xend=round(sim_row)-params_symb.min(params_symb.xoind);
                end
            else   
                if newrow
sim_row=(params_cont.Ad*xl+params_cont.Bd*ul+params_cont.minp);
round(sim_row);
params_symb.min(params_symb.xoind);
xend=round(sim_row)-params_symb.min(params_symb.xoind);
                    newrow=0;
                else
sim_row=sim_row+params_cont.Ad(:,1);
round(sim_row);
params_symb.min(params_symb.xoind);                    
xend=round(sim_row)-params_symb.min(params_symb.xoind);
                end
            end
            
            %Out-of-bounds Check
cond1=sum((xend(:)>params_symb.num(params_symb.xoind))+(xend(:)<zeros(params_symb.n,1)))==0;
            if(cond1) %If the transition is within bounds add it as is
trans(:,i)=[xuo; xend];
            else
trans(:,i)=[xuo; xoob]; 
            end
trans(:,i);
        else
        %% NON-DETERMINISTIC LINEAR SYSTEMS

            if isfield(params_cont,'num_sw')
                if newrow
                    lmaxtemp_row=params_cont.Ad{ul(end)+1}*xl+params_cont.Bd{ul(end)+1}*ul(1:end-1)+params_cont.minp{ul(end)+1}+params_cont.vmax{ul(end)+1};
                    lmintemp_row=params_cont.Ad{ul(end)+1}*xl+params_cont.Bd{ul(end)+1}*ul(1:end-1)+params_cont.minp{ul(end)+1}+params_cont.vmin{ul(end)+1};
                    newrow=0;
                else
                    lmaxtemp_row=lmaxtemp_row+params_cont.Ad{ul(end)+1}(:,1);
                    lmintemp_row=lmintemp_row+params_cont.Ad{ul(end)+1}(:,1);
                end

            else
                
                if newrow
lmaxtemp_row=params_cont.Ad*xl+params_cont.Bd*ul+params_cont.minp+params_cont.vmax;
lmintemp_row=params_cont.Ad*xl+params_cont.Bd*ul+params_cont.minp+params_cont.vmin;
                    newrow=0;
                else
lmaxtemp_row=lmaxtemp_row+params_cont.Ad(:,1);
lmintemp_row=lmintemp_row+params_cont.Ad(:,1);
                end

            end
            % This is the check to see if the intersection with neighbouring cells
            % is empty or not.
            for k=1:params_symb.n
                if(mod(lmaxtemp_row(k),1)~=1/2)
                    lmaxtemp(k)=round(lmaxtemp_row(k));
                else
                    lmaxtemp(k)=floor(lmaxtemp_row(k));
                end;
                if(mod(lmintemp_row(k),1)~=1/2)
                    lmintemp(k)=round(lmintemp_row(k));
                else
                    lmintemp(k)=ceil(lmintemp_row(k));
                end;
            end;
lmaxtemp;
lmintemp;
			% Shift to get labels between 0 and maximum label
lmax=lmaxtemp-params_symb.min(params_symb.xoind);
lmin=lmintemp-params_symb.min(params_symb.xoind);
            %%%

            %%% Add all the transitions generated:
            x=lmin;
            ntrans=1;
            jl=1;

            % The 3 lines below are a fix needed in the very special case
            % in which both lmaxtemp and lmintemp line in the boundary of
            % two cells. The line below resolves the ambiguity of the end
            % point to take. In any case, this should not happen unless
            % the matrix A is very ill-conditioned.
            if lmin(params_symb.n)>lmax(params_symb.n)
                lmax=lmin;
            end;
            % End-fix %
            
            while(jl<=params_symb.n && x(params_symb.n)<=lmax(params_symb.n))  
            jl=1;
                xendtemp(:,ntrans) = x;
                while(jl<params_symb.n && x(jl)>=lmax(jl))
                    x(jl)=lmin(jl);
                    jl=jl+1;
                end
                x(jl)=x(jl)+1;
                ntrans=ntrans+1;
            end
            
xend=xendtemp(:,1:ntrans-1);
        
            %Add transitions (non-deterministic) all in one shot
            if ~isempty(xend) % xend should never be empty... but if it is, we don't include that transition
                numtrans=size(xend,2);
            %Out-of-bounds Check
cond1=sum((max(xend,[],2)>params_symb.num(params_symb.xoind))+(min(xend,[],2)<zeros(params_symb.n,1)))==0;
                if(cond1) %only if none of the transitions went out of bounds
newtrans=[kron(xuo,ones(1,numtrans)); xend];
                    trans_tmp(:,counter:counter+numtrans-1)=newtrans;        
                    counter=counter+numtrans;
                else
newtrans=[xuo; xoob];
					trans_tmp(:,counter)=newtrans;
                    counter=counter+1;
                end
            end;
        end;
        
    else %NON-LINEAR SYSTEMS
        xend=pss_dynamics(params_symb,params_cont,xuo(params_symb.xoind),xuo(params_symb.uind));

        %Add transitions (non-deterministic) all in one shot
        if ~isempty(xend) % xend should never be empty... but if it is, we don't include that transition
        %Out-of-bounds Check
            cond1=sum((max(xend,[],2)>params_symb.num(params_symb.xoind))+(min(xend,[],2)<zeros(params_symb.n,1)))==0;
            if(cond1) %only if none of the transitions went out of bounds
                newtrans=[kron(xuo,ones(1,size(xend,2))); xend];
                trans=[trans newtrans];        
            else
                trans=[trans [xuo; xoob]]; 
            end
        end;
    end;

    % Scan through labels
    %Increase the labels starting from the leftmost in the vector [xini u]
    j=1;
    while (j<=npm && xuo(j)>=params_symb.num(j))
        xuo(j)=0;
        j=j+1;
    end
    if j>npm
            break;
    else
        xuo(j)=xuo(j)+1;
    end
    if j~=1
        newrow=1;
    end;
    
end

%% RETURN

if ~isempty(trans_tmp) && isfield(params_cont,'Ad') && isfield(params_cont,'Bd') && params_symb.deter==2 && ~params_cont.custom %NON-DETERMINISTIC LINEAR SYSTEMS NOT CUSTOM
    trans=trans_tmp(:,1:counter-1);
end;

xnext=xuo(params_symb.xoind);
unext=xuo(params_symb.uind);
