% 1. note that the grid is spanned only within the
%    [-1 1]^n. Thus, everything needs to be scaled, such 
%    that the whole state space is contained within [-1 1]^n

% 2. the polytopes are assumed to be given in form of
%    
%     P={ x | Hx <= h }
% 
%     where the matrix H and vector h are provided


%
% In order to run the exmaple,
% you need to have the MPT toolbox in your MATLAB path
%

close all; 
clear set;



% Example of a 2D system with one input
% and initial set is given by two triangles

% first polytope
H1 = [-0.9304   0.3666;
      -0.0811  -0.9967;
       0.6835   0.7299];
h1 = [ 0.0940; -0.3; .71]-H1*[0.7; 0.7];


% second polytope
H2 = [ 0.8835    0.4684;
      -0.8742   -0.4855;
      -0.6475    0.7621;
       0.0608   -0.9982;
       0.3634   -0.9316];
h2 = [ 0.9761; -0.5748; 0.2115; -0.2590; 0.0091];

%H=[-1 0 0;
%    1 0 0;
%    0 1 0;
%    0 -1 0;
%    0  0 1;
%    0  0 -1];
%
%h=[2/10 2/10 6/10 6/10 4/10 4/10]';

% mpt toolbox at hand? plot polytopes...
if 1
  hold on
  plot(polytope(H1,h1),'g')
  plot(polytope(H2,h2),'g')
end

% state space
dim=2;

% init cudd package 
set=mcis(dim);

% add first polytope to bdd
set.addPolytope(H1,h1);

% add second polytope to bdd
set.addPolytope(H2,h2);

% compute mcis with 8 bits resolution in each dimension
tic
set.compute(8);
toc


% get nodes of the outer approxiamtion
nodes=set.getNodes(1);
plot(nodes(:,1),nodes(:,2),'k.')
hold on
axis([-1 1 -1 1])

%% get nodes of the inner approxiamtion
nodes=set.getNodes(0);
plot(nodes(:,1),nodes(:,2),'bo')
axis([-1 1 -1 1])
