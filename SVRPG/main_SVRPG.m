%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On the Oracle Complexity Reduction of the Linear Quadratic Regulator 
%    (LQR) via Stochastic Variance-Reduced Policy Gradient (SVRPG)
%              Leonardo F. Toso, Han Wang, James Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear all; close all

load('workspace.mat')
load("U2.mat")

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


%% LQR via PG with two-point gradient estimation


%Inner loop - Smoothing radius
r=5e-2;

%Number of outer-loop rollouts

ns=50;

%Number of inner-loop rollouts

ns_hat=25;


%Number of realizations for the confidence bounds (95% confidence)

nr=5;

%Step-size:

eta=1e-4;


%Number of epochs

N=125;

%Epoch length

T=4;


%Number of set of matrices 

q=1;

costs_SVRPG=zeros(nr,N*T);
grad_norm_SVRPG=zeros(nr,N*T-1);
for i=1:nr
    [costs_SVRPG(i,:),grad_norm_SVRPG(i,:)]=LQR_SVRPG(A,B,Q,R,K0,N,T,r,ns,ns_hat,eta,initial_gap,cost_optimal,cost_initial,q,U,U2);
end

