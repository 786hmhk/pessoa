
classdef mcis_mpt < handle

  methods 

    %% constructor 
    function this=mcis_mpt(A,B,D,H)
      
      % discrete state vector
      Q=cell(size(D,1),1);

      this.Q=Q;
      this.A=A;
      this.B=B;
      this.D=D;

      this.Pu=[];

    end

    function this=addOutput(i,P)
      
      this.Q{i}=P

    end


    function this=addInput(Pu)
      
      this.Pu=Pu;

    end

    %%
    function display(this)
    % display function
      disp(' ')
      disp('maximal controlled invariant set via polyhedral approach');
      disp(' ')
    end

    %% compute invariant set
    function compute(this)

      %% Prepare iteration
      opt.abs_tol=1e-9;
      opt.rel_tol=1e-9;
      opt.verbose=2;
      opt.facecolor=[0.2 0.4 0.2];
      opt.projection=6;

      A=this.A;
      B=this.B;


      Q=this.Q;
      Pu=this.Pu;

      n=length(A);
      Nq=length(Q);
      Nm=length(Pu);


      % start
      iter=0;

      bb=0;

      while iter<=100

        marker=0;
        maker2=0;
        for i=1:Nq

          P=Q{i};

          if isempty(P) || ~isfulldim(P)
            bb=bb+1;
            continue;
          end

          [H h]=double(P);

          if Nm
            [Hu hu]=double(Pu);
            if Nm==1,
              Hu={Hu}; 
              hu={hu};
            end
          end

          if length(P)==1,
            H={H}; 
            h={h}; 
          end

          clear PP;
          if Nm
            for i=1:length(P)
              for j=1:Nm
                Ht=blkdiag(H{i},Hu{j});
                ht=[h{i};hu{j}];

                PP((i-1)*Nm+j)=projection(polytope(Ht*[A B; zeros(m,n) eye(m)],ht),1:n,opt);
              end
            end
          else
            for i=1:length(P)
                Ht=H{i};
                ht=h{i};

                PP(j)=projection(polytope(Ht*[A B],ht),1:n,opt);
              end
            end
          end
           
          for j=1:Nq
            if D(i,j)==1
              PP=polytope(intersect(Q{j},PP,opt)); %CPre(P) * P

              if isempty(P) || ~isfulldim(P)
                Q{j}=[];
              else
                Q{j}=PP;
              end
            
            end
          end

            if ~le(P,PP,opt)
              marker=1;

           %   if plotting
           %     plot(projection(PP,[1 2 3]),opt);
           %     drawnow
           %   end
               break;
            end 
             
          end
        end

        if marker==0
          disp(['P_i \subset P_i+1'])
          break;
        end 
          
        if bb==Nq
          display('The resulting polytope is empty');
          break;
        end

        iter=iter+1;
        disp(['iteration ',num2str(iter)])
      end

    end


  end % end methods

  properties (GetAccess = 'public', SetAccess='private')

    % system data
    A;
    B;
    % discrete data
    Q; % states
    D; % transitions
    H; % output function
    
    % atomic propositions
    P;

    % input constraints
    Pu;
  end % end private properties



end
