
classdef mcis < handle

  methods 

    %% constructor 
    function this=mcis(dim,varargin)
    %
    % Approximated computation of the maximal control invariant set 
    % 
    % The state space is fixed to [-1 1]^n. All sets need 
    % to be scaled such that they fit into [-1 1]^n.
    % 
    % The system is assumed to be in brunovsky normal form
    % 
    % SET=MCIS(DIM,VARARGIN)
    %
    % dim   - state space dimension
    % 
    % varargin
    % 'mu'  - vector of ctrb indices (default mu=n)
    % 'addRoundingPrecision' - rounding precision in the add leafs
    %

      % if there exist another mcis class
      % leave libraries loaded
      s=evalin('base','whos');
      if length(strfind([s(:).class],'mcis'))
        error(['Only ONE instance of a mcis class supported!'...
                 ' Clear all other mcis instances first AND unload mcis library'])
      end

      if libisloaded('mcis')
        try 
          unloadlibrary('mcis')
        catch
          error(['Library mcis is still loaded'...
                 'Please manually unload mcis library'])
        end
      end

      if (nargin<1)
        error('Input arguments wrong: Please see help mcis')
      end

      dim=int32(dim);

      % ctrb indices default value (one input)
      mu=dim;

      pre=1e-12;

      % process input arguments
      if ~isempty(varargin)
        for i=1:2:nargin-2
          switch lower(varargin{i})
          case 'mu'
            mu=int32(varargin{i+1});
          case 'addroundingprecision'
            pre=double(varargin{i+1});
          otherwise
            error('Input arguments wrong: Please see help mcis')
          end
        end
      end

      if (sum(mu)~=dim)
        error('Controllability indices must sum up to dim')
      end


      % default suffix
      switch computer('arch')
        case{'glnx86','glnxa64'}
          suffix='.so';
        case{'maci','maci64'}
          suffix='.dylib';
        otherwise
          error('System Architecture not supportet.')
      end
   
      % load c library (like mex files)
      warning off MATLAB:loadlibrary:TypeNotFoundForStructure
      [~, ~]=loadlibrary(['mcis',suffix],'mcis.h');

      % inputs (ctrbl indices)
      Nu=int32(length(mu));
      %mu=libpointer('int32Ptr',mu);

      % fill system data
      sys.Nu=Nu;
      sys.mu=mu;
      sys.dim=dim;
      sys.Nbits=1;
      % leaf node values are rounded 
      sys.precision=pre; 

      % create pointer
      sys_ptr=libpointer('sys_Ptr',sys);

      % initialize decision diagrams
      calllib('mcis','mcis_init',sys_ptr);
      % store things in class
      this.suffix=suffix;
      this.sys=sys_ptr;
      this.intptr=libpointer('int32Ptr',zeros(dim,1));

      % grid resolution
      this.dx=0;

    end

    %% deconstructor
    function delete(this)

      if libisloaded('mcis')
        refcount=calllib('mcis','mcis_close');
        if refcount
          error('Reference count should be zero')
        end
        
        % free pointer
        %mu=this.sys.mu; this.sys.mu=[]; clear mu;
        sys=this.sys; 
        this.sys=[]; 
        clear sys;

        intptr=this.intptr; 
        this.intptr=[];
        clear intptr;

        unloadlibrary('mcis')
      end
    end

    %%
    function display(this)
    % display function
      disp(' ')
      disp('maximal controlled invariant set data structure');
      disp(' ')
    end

    %% compute invariant set
    function compute(this,Nbits)
    % Compute the invariant set 

      calllib('mcis','mcis_compute',int32(Nbits));

      % update grid resolution
      this.dx=2^(-Nbits+1);

    end

    %% compute bdd with points inside the polytope
    %
    function addPolytope(this,H,h)
    %
    % this.addPolytope(H)
    %
    % updates the bdd such that it containts
    % the grid points inside P={ x | Hx <= 1} \subseteq B(0,1)
    % 
    % grid points spann [-1 1]^n
    % with 2^Nbits points in each dimension
    %

      
      res=calllib('mcis','mcis_addPolytope',H,h,int32(size(H,1)));

      if (res==0)
        error('Adding polytope failed!')
      end

    end

    %% 
    function x=query(this,x)
    % C=THIS.QUERY(x) finds hypercube containing x
    % C is empty if x is not contained in the set
    % 


      x=x(:)';
      if any(size(x))~=1
        error('Please provide correct dimensions.')
      end

      if (max(x)>1) || (min(x)<-1)
        error('Query point must be within the unit cube.')
      end

      dim=this.sys.v.dim;
      Nbits=this.sys.v.Nbits;

      y=libpointer('int32Ptr',floor((1+x).*2^(Nbits-1)));

      res=calllib('mcis','mcis_query',y);
      
      % successful query
      if res 
        y=double(y.v);
        x=y/2^(Nbits-1)-1;
      else
        x=[];
      end


    end

    % count cubes in bdd
    function N=getNumCubes(this,io)

      % number of nodes (up to don't cares)
      N=calllib('mcis','mcis_getNumCubes',io);

    end
 
    % plot nodes inside the set
    function nodes=getNodes(this,io,varargin)
    % GETNODES(THIS,IO) returns the nodes that are 
    % contained in the set, where the set corresponds to
    % the inner (IO=0) or outer (IO=1) approximation
    % 
    % GETNODES(THIS,IO,PRO) returns the nodes
    % contained in the set projected onto the dimensions
    % provide by the vector PRO
    % 

      dim=this.sys.v.dim;
      Nbits=this.sys.v.Nbits;

      m=1;


      % number of nodes (up to don't cares)
      N=calllib('mcis','mcis_getNumCubes',io);

      if (N==0)
        nodes=[];
        return;
      end


      N=ceil(N/m);

      % nodes encoded as binary numbers
      % each row corresponds to one node
      bddCubes=repmat('0',dim*Nbits*N,1);
      bddCubes=calllib('mcis','mcis_getNodes',bddCubes,io,m);
      bddCubes=reshape(bddCubes,dim*Nbits,N)';

      % change don't care bits, indicated by '2' 
      % in the cube to 0 and 1
      for i=1:dim*Nbits
        idx=find(bddCubes(:,i)=='2');
        if ~isempty(idx)
          tmp=bddCubes(idx,:);
          tmp(:,i)='0';

          bddCubes(idx,i)='1';
          bddCubes=[bddCubes;tmp];
        end
      end

      % determine the dimensions 
      % for which the nodes are computed
      if isempty(varargin)
        pro=1:dim;
        n=dim;
      else
        pro=varargin{1};
        n=length(pro);
      end


      % convert binary encoding to dec
      nodes=zeros(size(bddCubes,1),n);
      for i=1:n
        j=pro(i);
        %nodes(:,i)=bin2dec(bddCubes(:,(j-1)*Nbits+1:j*Nbits))/2^Nbits;
        nodes(:,i)=bin2dec(bddCubes(:,j:dim:dim*Nbits))/2^(Nbits-1)-1;
      end

     
    end

    % compute upper bound on error in 
    % terms of the volume 
    function e=getError(this,varargin)
    % GETERROR(THIS,BOUND) computes the abusolute difference
    % between the volume of the inner and outer 
    % approximation 
    % 
    % for BOUND == 1 we get an upper bound on the error
    % 
    % error <= ( vol(OUTER) - vol(INNER) )/ vol(INNER)
    %
    % for BOUND == 0 we get a lower bound on the error
    %
    % error <= ( vol(OUTER) - vol(INNER) )/ vol(OUTER)
    %
    if isempty(varargin)
      error('Please chose upper or lower bound.')
      error('see help this.getError.')
    end

      e=calllib('mcis','mcis_getError',varargin{1})*100;

    end
 
    % some information about the dd 
    function getInfo(this,varargin)
    % INFO(THIS,P) prints some information about the 
    % becision diagrams. If P==1, a dot file containing
    % the approximations is written.

      % print the dot file or not
      if isempty(varargin)
        pr=0;
      else
        pr=1;
      end
 
      % print info
      calllib('mcis','mcis_getInfo',pr);

    end
 

  end % end methods

  properties (GetAccess = 'public', SetAccess='private')

    % suffix for different arch (.dylib, .so)
    suffix;
    % system data
    sys;
    % allocated pointer (for permformance reasons)
    intptr;

    % grid resolution
    dx;

  end % end private properties



end
