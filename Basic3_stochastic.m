%%      AUTHORS:    Mirka Mandich, Seattle University
%                   Tianwei Xia, University of Tennessee, Knoxville
%
%       MENTORS:    Dr. Kai Sun, University of Tennessee, Knoxville
%                   Dr. Kevin Tomsovic, University of Tennessee, Knoxville
%
%       PROJECT:    Optimal PMU Placement Based on Robust Optimization 
%                   CURENT REU 2018
%
%   DESCRIPTION:    Optimizes locations of PMUs in IEEE 9-Bus System 
%                   (considers zero-injection buses)
%                   (cost-optimized; lamda = 0)
%                   (uses stochastic techniques)
%
%    REFERENCES:    [4] Optimal PMU Placement based on Integer Linear Programming
%                   [3] Contingency-constrained PMU Placement in Power
%                   Networks
%                   [5] Chance-constrained
%
%       UPDATED:    7/19/2018, 3:00 PM
%%        INPUT:
clc
clear
database = [];

%confidence interval
disp('What confidence interval for PMU success would you like? (For example: "0.9" ensures 90%.)');
user = input('Enter a value from 0 to 1: ');
if isempty(user)
    user = 0.95;
    disp('Default value used: 0.95');
    disp(' ');
end

%iterations
disp('How many tests would you like to run? (For example: "100".)');
user2 = input('Enter a greater than 10: ');
if isempty(user2)
    user2 = 50;
    disp('Default value used: 50');
    disp(' ');
end

%%    STOCHASTIC:
%Eta_Failure contains probability of branch failures from [7], Table I
%0 indicates no connection or "self-connected" PMUs (diagonal)
Eta_Failure = zeros(9);
Eta_Failure(1,4) = 0.0089;
Eta_Failure(2,7) = 0.0033;
Eta_Failure(3,9) = 0.0055;
Eta_Failure(4,1) = 0.0089; Eta_Failure(4,5) = 0.0054; Eta_Failure(4,6) = 0.0046;
Eta_Failure(5,4) = 0.0054; Eta_Failure(5,7) = 0.0044;
Eta_Failure(6,4) = 0.0046; Eta_Failure(6,9) = 0.0091;
Eta_Failure(7,2) = 0.0033; Eta_Failure(7,5) = 0.0044; Eta_Failure(7,8) = 0.0053;
Eta_Failure(8,7) = 0.0053; Eta_Failure(8,9) = 0.0023;
Eta_Failure(9,3) = 0.0055; Eta_Failure(9,6) = 0.0091; Eta_Failure(9,8) = 0.0023;

%Eta_Failure_Tri removes redundancies (internal use only)
Eta_Failure_Tri = tril(Eta_Failure);

%Eta_Success
[row, col] = size(Eta_Failure);
Eta_Success = zeros(row);
for r = 1:row
    for c = 1:col
        if Eta_Failure(r, c) ~= 0
            Eta_Success(r, c) = 1-Eta_Failure(r, c);
        end
    end
end

%Eta_Success_Tri
Eta_Success_Tri = tril(Eta_Success);

[row, col] = find(Eta_Failure_Tri); %finds locations of non-zeros in Eta (AKA the bus connections)
for i=1:length(row) %for every connection...
            %figure();
            %title(['gaussian distribution data (' num2str(i) ') for branch ' num2str(row(i)) ' to ' num2str(col(i))]);
            %hold on;
    
    mu = Eta_Failure_Tri(row(i),col(i)); %center of gaussian plot; mean
    x = [0:.005:1]; %range that plots every 0.005
    sigma = 0.05; %standard deviation
    r = abs(normrnd(mu,sigma,[1,user2]));
    
    for j=1:10 %...create monte carlo values
        norm = normpdf(x,r(j),sigma);
            %plot(x,norm)
    end
            %hold off;

    %Eta_Failure_Gaussian
    Eta_Failure_Gaussian(i,:) = r; 
end

%Eta_Success_Gaussian
[row2, col2] = size(Eta_Failure_Gaussian); %col2 = branch connection (1:9); row2 = monte carlo test (1-user2)
Eta_Success_Gaussian = zeros(row2, col2);
for r2 = 1:row2
    for c2 = 1:col2
            Eta_Success_Gaussian(r2, c2) = 1-Eta_Failure_Gaussian(r2, c2);
    end
end

%%  OPTIMIZATION:    (ONCE) WITHOUT MONTE CARLO TESTING
%inequality constraints (A*x =< b)
A = -Eta_Success+-eye(size(Eta_Success)); %[4]eq4 x-values
A2 = zeros(9,8); %[4]eq4 y-values
A2(2,1) = A(2,7); %y27 %z7
A2(5,2) = A(5,7); %y57
A2(8,3) = A(8,7); %y78
A2(7,4) = A(7,7); %y77
A2(3,5) = A(3,9); %y39 %z9
A2(6,6) = A(6,9); %y69
A2(8,7) = A(8,9); %y89
A2(9,8) = A(9,9); %y99
A = [A A2];
b = -(ones(1,9))'*user; %our confidence interval

%equality constraints (Aeq*x = beq); zero-injection
Aeq = []; 
Aeq = zeros(2,9);
Aeq2 = [1 1 1 1 0 0 0 0; %[4]eq4 z7
        0 0 0 0 1 1 1 1]; %[4]eq4 z9
% test1 = nonzeros(A2);
% Aeq2 = [-test1(1:4).' 0 0 0 0; %[4]eq4 z7
%         0 0 0 0 -test1(5:8).']; %[4]eq4 z9
Aeq = [Aeq Aeq2]; %zero injection
beq = [1 1]; %[4]eq4 z7 and z9 must equal 1
%beq = [user2 user2]; %[4]eq4 z7 and z9 must equal 1

f = [ones(1,9) zeros(1,8)]; %In Basic1, optimized f for cost and TMR. Here, only cost matters (lamda=0)
intcon = 1:17; %all must be integers
lb = zeros(1,17); %enforces binary
ub = ones(1,17); %enforces binary

[x,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
disp(x) %displays a PMU setup which maintains full observability with optimal cost
disp(sum(x)); %number of PMUs in setup

%%  OPTIMIZATION:   WITH MONTE CARLO TESTING
[row3, col3] = size(Eta_Success_Gaussian); %(9, user2)
[row4, col4] = size(Eta_Success_Tri); %(9, 9)
pmu = [];

for c3 = 1:col3 %run it user2 times
    r3 = 1;
    for r4 = 1:row4
        for c4 = 1:col4
            if Eta_Success_Tri(r4, c4) ~= 0 %if there's a connection...
                Eta_Success_Tri(r4, c4) = Eta_Success_Gaussian(r3, c3); %...replace with monte carlo...
                r3 = r3+1; %...and move onto the next one
            end
        end
    end
    Eta_Success_Tri;
    Eta_Success = Eta_Success_Tri+rot90(flip(Eta_Success_Tri),-1);
    A = -Eta_Success+-eye(size(Eta_Success));
    A2 = zeros(9,8); %[4]eq4 y-values
    A2(2,1) = A(2,7); %y27 %z7
    A2(5,2) = A(5,7); %y57
    A2(8,3) = A(8,7); %y78
    A2(7,4) = A(7,7); %y77
    A2(3,5) = A(3,9); %y39 %z9
    A2(6,6) = A(6,9); %y69
    A2(8,7) = A(8,9); %y89
    A2(9,8) = A(9,9); %y99
    A = [A A2]; %our new A matrix!

    [x,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
    disp(x); %displays a PMU setup which maintains full observability with optimal cost
    pmu = [pmu, sum(x)];
end
disp(pmu); %number of PMUs for each setup

%%          PLOT:
figure();
bar(1:user2, pmu, 'FaceColor', [0.8 0.9 0.8], 'EdgeColor', [0.5 0.9 0.5]);
xlabel(['Results from ' num2str(user2) ' Monte Carlo Simulations']);
ylabel('Number of PMUs Required per Simulation');
grid on;
title('Cost-Optimized Results for Model 2.1');
hold on; plot(xlim,[mean(pmu) mean(pmu)], 'r');
%txt1 = ['The average number of PMUs needed will be ' num2str(mean(pmu))];
text(1,0.5,['PMUs needed, on average = ' num2str(mean(pmu))]);

disp(['PMUs needed, on average = ' num2str(mean(pmu))]);
