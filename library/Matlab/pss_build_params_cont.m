function params_cont=pss_build_params_cont(params_symb,A,B,C,dset)

%
% params_cont=pss_build_params_cont(params_symb,A,B,C,dset)
%
% Constructor for the data structure "params_cont" containing the continuous-time
% model information necessary to construct the abstraction. The dynamics of the 
% continuous model is assumed to be of the form $\dot{x}=Ax+Bu+Cd$, where u is a 
% control input and d a disturbance input
%   
% INPUTS: params_symb - structure containing information about the abstraction
%                       to be constructed;
%         A           - See above;
%         B           - See above;
%         C           - See above;
%         dset        - matrix defining the uncontrolled input set, each row 
%                       contains the minimum and maximum value for that coordinate.
%
% OUTPUTS: params_cont - resulting structure.
%
% For more details see the online documentation.
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab 2009.

% params_cont.Ad
% params_cont.Bd
% params_cont.minp
% params_cont.vmax
% params_cont.vmin
% params_cont.Bnd
% params_cont.vmin
% params_cont.vmax
if nargin >=3
     if isfield(params_symb,'num_sw')
         params_cont.num_sw=params_symb.num_sw;
         if nargin>3
%            for i=1:params_cont.num_sw
%             [foo,nb] = size([B{i},C{i}]);
%             s = expm([[A{i} B{i}]*params_symb.tau; zeros(params_symb.m,params_symb.n+params_symb.m)]);
%             params_cont.Bnd{i}=s(1:params_symb.n,params_symb.n+params_symb.m+1:params_symb.n+nb)/params_symb.eta; 
%             params_cont.Ad{i} = s(1:params_symb.n,1:params_symb.n);
%             params_cont.Bd{i} = s(1:params_symb.n,params_symb.n+1:params_symb.n+params_symb.m)*(params_symb.mu/params_symb.eta);
%            end
            error('Switched systems with disturbances not yet supported.')
         else
           for i=1:params_cont.num_sw
            [foo,n] = size(A{i});
            [foo,nb] = size(B{i});
            s = expm([[A{i} B{i}]*params_symb.tau; zeros(nb,n+nb)]);
            params_cont.Ad{i} = s(1:n,1:n);
            params_cont.Bd{i} = s(1:n,n+1:n+nb)*(params_symb.mu/params_symb.eta);
           end
        end
    else
        if nargin>3
            [foo,nb] = size([B,C]);
            s = expm([[A B C]*params_symb.tau; zeros(nb,params_symb.n+nb)]);
            params_cont.Bnd=s(1:params_symb.n,params_symb.n+params_symb.m+1:params_symb.n+nb)/params_symb.eta; 
            params_cont.Ad = s(1:params_symb.n,1:params_symb.n);
            params_cont.Bd = s(1:params_symb.n,params_symb.n+1:params_symb.n+params_symb.m)*(params_symb.mu/params_symb.eta);
        else
            s = expm([[A B]*params_symb.tau; zeros(params_symb.m,params_symb.n+params_symb.m)]);
            params_cont.Ad = s(1:params_symb.n,1:params_symb.n);
            params_cont.Bd = s(1:params_symb.n,params_symb.n+1:params_symb.n+params_symb.m)*(params_symb.mu/params_symb.eta);
        end
    end
else
    error('insuficient input parameters. See help.');
end

if isfield(params_symb,'num_sw')
    min_all_inputs=params_symb.min(params_symb.uind);
    min_only_inputs=min_all_inputs(1:end-1);
    
    for i=1:params_symb.num_sw
        params_cont.minp{i}=params_cont.Ad{i}*params_symb.min(params_symb.xoind)+params_cont.Bd{i}*min_only_inputs;
        params_cont.vmax{i}=sum(abs(params_cont.Ad{i}),2)*0.5;
        params_cont.vmin{i}=-params_cont.vmax{i};
    end
else
     
    params_cont.minp=params_cont.Ad*params_symb.min(params_symb.xoind)+params_cont.Bd*params_symb.min(params_symb.uind);
    params_cont.vmax=sum(abs(params_cont.Ad),2)*0.5;
    params_cont.vmin=-params_cont.vmax;

    if nargin>3
        dset_center=0.5*(dset(:,1)+dset(:,2));
        dset_var=dset(:,2)-dset_center;
        params_cont.vmin=params_cont.Bnd*dset_center+params_cont.vmin-sum(abs(params_cont.Bnd*diag(dset_var)),2); % Missing effect of input u_d not matching u_c (gamma(\mu)).
        params_cont.vmax=params_cont.Bnd*dset_center+params_cont.vmax+sum(abs(params_cont.Bnd*diag(dset_var)),2);
    end

end