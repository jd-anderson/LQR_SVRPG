%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On the Oracle Complexity Reduction of the Linear Quadratic Regulator 
%    (LQR) via Stochastic Variance-Reduced Policy Gradient (SVRPG)
%              Leonardo F. Toso, Han Wang, James Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear all; close all


%% System matrices 

A=[1.20 0.50 0.40;
   0.01 0.75 0.30;
   0.10 0.02 1.50];

nx=size(A,1);
nu=1;

B=[0.25;1;1/2];


%% Initial stabilizing feedback controller gain 

K0=[0.15   -0.45    3.80];

%Cost matrices

Q=2*eye(nx);
R=(1/2)*eye(nu);

cost_initial=cost(A,B,Q,R,K0);


%% Optimal controller and optimal cost

K_opt=dlqr(A,B,Q,R);

cost_optimal=cost(A,B,Q,R,K_opt);

%% Initial gap

initial_gap=cost_initial-cost_optimal;



%% LQR via PG with known model (model-based)


% Number of iterations

N=500;


% Step-size

eta=1.5e-4;


% Number of realizations

nr=1;

costs_model_based=zeros(nr,N);
grad_norm_model_based=zeros(nr,N-1);
for i=1:nr
    [costs_model_based(i,:),grad_norm_model_based(i,:)]=LQR_model_based(A,B,Q,R,K0,N,eta,initial_gap,cost_optimal,cost_initial);
end

save('costs_model_based.mat','costs_model_based');
save('grad_norm_model_based.mat','grad_norm_model_based');
