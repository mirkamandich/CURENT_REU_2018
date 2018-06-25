%       AUTHORS:    Mirka Mandich, Seattle University
%                   Tianwei Xia, University of Tennessee, Knoxville
%
%       MENTORS:    Dr. Kai Sun, University of Tennessee, Knoxville
%                   Dr. Kevin Tomsovic, University of Tennessee, Knoxville
%
%       PROJECT:    Optimal PMU Placement Based on Robust Optimization 
%                   Considering N-k Cascading Failures
%                   CURENT REU 2018
%
%   DESCRIPTION:    Optimizes locations of PMUs in IEEE 9-Bus System 
%                   (does not consider zero-injection buses)
%
%    REFERENCES:    [4] Optimal PMU Placement based on Integer Linear Programming
%                   [3] Contingency-constrained PMU Placement in Power
%                   Networks
%
%       UPDATED:    6/25/2018, 11:15 AM

clc
clear
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1261\cplex\matlab'); %CPLEX solver
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1261\cplex\matlab\x64_win64');
database=[];

f1=ones(1,9); %cost; want to minimize for fewest PMUs (each PMU "costs" 1)
f2=[2 2 2 4 3 3 4 3 4]; %redundancy; want to maximize (AKA number of connections to other buses, TMR [4]eq7)

for lamda=0:0.01:1 %weighting; 0 = only cost matters, 1 = only redundancy matters
f=f1*(1-lamda)-lamda*f2; %objective function; [4]eq9 where w = lambda, x_i = f1, TMR = f2

%inequality constraints (A*x =< b)
A=zeros(9,9); %[4]eq2 x-values
A(1,1)=-1;A(1,4)=-1;
A(2,2)=-1;A(2,7)=-1;
A(3,3)=-1;A(3,9)=-1;
A(4,4)=-1;A(4,1)=-1;A(4,5)=-1;A(4,6)=-1;
A(5,5)=-1;A(5,4)=-1;A(5,7)=-1;
A(6,6)=-1;A(6,4)=-1;A(6,9)=-1;
A(7,7)=-1;A(7,2)=-1;A(7,5)=-1;A(7,8)=-1;
A(8,8)=-1;A(8,7)=-1;A(8,9)=-1;
A(9,9)=-1;A(9,3)=-1;A(9,6)=-1;A(9,8)=-1;
b=-(ones(1,9))'; %no function may equal -1; full observability

%equality constraints (Aeq*x = beq); null because not considering zero-injection
Aeq=[];
beq=[];

lb=zeros(1,9); %enforces binary
ub=ones(1,9);  %enforces binary
intcon = 1: size(A); %all must be integers

%options = optimoptions(@intlinprog,'OutputFcn',@savemilpsolutions,'PlotFcn',@optimplotmilp);
%x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
[x,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
disp(x) %displays a PMU setup which maintains full observability with various cost or redundancy
Tempdata=[lamda x' f1*x f2*x fval];
database=[database;Tempdata];
end

database(1,12);
plot(database(:,1),database(:,11),'-.',database(:,1),database(:,12),'-',database(:,1),database(:,13),'.')
grid on;
xlabel('w')
ylabel('object value')
legend('Minimum number of PMU','Maximum number of redundancy','Combination result')
title('multi-objective programming for model 1.0')