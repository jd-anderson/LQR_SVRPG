%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On the Oracle Complexity Reduction of the Linear Quadratic Regulator 
%    (LQR) via Stochastic Variance-Reduced Policy Gradient (SVRPG)
%              Leonardo F. Toso, Han Wang, James Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear all; close all

load('workspace.mat')
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


%% LQR via PG with ZO and two point gradient estimation 


%Smoothing radius

r=1e-4;

%Number of trajectories for the full gradient estimation

ns=13;

%Number of realizations for the confidence bounds (95% confidence)

nr=5;

%Number of iterations

N=500;

%Step-size

eta=1e-4;


costs_ZO2P=zeros(nr,N);
grad_norm_ZO2P=zeros(nr,N-1);

for i=1:nr
    [costs_ZO2P(i,:),grad_norm_ZO2P(i,:)]=LQR_PG_ZO2P(A,B,Q,R,K0,N,r,ns,eta,initial_gap,cost_optimal,cost_initial,U);
end

save('costs_ZO2P.mat','costs_ZO2P');
save('grad_norm_ZO2P.mat','grad_norm_ZO2P');
