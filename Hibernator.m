clear all;
clc;
close all;

addpath(' C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64')
% model for hibernator
%% Part I) parameters for input 
n = 10; % disks
J = 3; % speeds (Off - Low Speed - High Speed) 
Rlimit = 1000 * 0.001; % msec response time limit
Tepoch = 30 * 60 ; % min epoch length

% The speed of the disks in the previous run, 10 Random integers between 1 and 10
% 1 Off 
% 2 Low Speed
% 3 High speed
Initial = randi([1 3],1,10);

PIdle(1:n,1) = 1; % Power of off when Idle
PIdle(1:n,2) = 3.5; % Power of low speed when Idle
PIdle(1:n,3) = 6.9; % Power of high speed when Idle

PActive(1:n,1) = 1; % Power of off when active
PActive(1:n,2) = 8.5; % Power of low speed when active
PActive(1:n,3) = 11.6; % Power of high speed when active

PTransition(1:n,1) = 1; % Power of transition to off
PTransition(1:n,2) = 4.7; % Power of transition to low speed
PTransition(1:n,3) = 8.9; % Power of transition to high speed

alpha = rand(1,n) * 3 + 1; % arrival rate for disk i - Randomely assign between 9 to 11 continous 

tij_exp(1:n,1) = 1000 ; % sec expected service time when off. Was set very large to give a high penalty
tij_exp(1:n,2) = 2 * 0.001; % msec expected service time low speed
tij_exp(1:n,3) = 1.4 * 0.001; % msec xpected service time high speed

tij_var(1:n,1) = 0; % sec variance of service time when off
tij_var(1:n,2) = 0.3; % sec variance of service time Low speed
tij_var(1:n,3) = 0.5; % sec variance of service time high speed

T(1:n,3) = 15; % sec transition period to high speed
T(1:n,2) = 11; % sec transition period to low speed
T(1:n,1) = 0; % sec transition period from any higher mode to any lower mode

Ni = ones(n); % number of request at disk i

%% Part II) Determining the Optimization parameters
rho = ones(n,J); % disk utilization
R1 = ones(n,J);
R2 = ones(n,J);
R = ones(n,J);
E = ones(n,J);
N = ones(1,n);

for i = 1:n
    for j = 1:J
       ExcInd = Initial(1,i); 
       % Average Response time 
       rho(i,j) = alpha(i)*tij_exp(i,j);
       N(i) = alpha(i)* Tepoch;
       R1(i,j) = T(i,j)*(1+rho(i,j))/2;
       R1(i,ExcInd) = 0; 
       R2(i,j) = alpha(i)*tij_exp(i,j) + alpha(i)^2*(tij_exp(i,j)^2+tij_var(i,j))/(2*(1-alpha(i)*tij_exp(i,j)));
       R(i,j) = (alpha(i)*R1(i,j)*T(i,j) + alpha(i)*(Tepoch-T(i,j))*R2(i,j))/(alpha(i)*Tepoch);
       R(i,ExcInd) = R2(i,ExcInd);
       % Energy 
       E(i,j) = PActive(i,j)*Tepoch*rho(i,j) + PIdle(i,j)*(Tepoch*(1-rho(i,j))-T(i,j)) + PTransition(i,j)*T(i,j);
       E(i,1) = PIdle(i,1)*Tepoch; 
       E(i,ExcInd) = PActive(i,ExcInd)*Tepoch*rho(i,ExcInd) + PIdle(i,ExcInd)*Tepoch*(1-rho(i,ExcInd));
    end    
end
TotalN = sum(N);
TotalReq = [N' N' N'] .* R;

%% Part III) The Cplex model

try
   Energy = reshape (E', n*J, 1);
   
   % Build model
   cplex = Cplex('Hibernator');
   
   % The objective is to minimize the energy consumption.
   cplex.Model.sense = 'minimize';
   
   obj   = Energy;
   lb    = zeros (n*J, 1);
   ub    = ones (n*J, 1);
   ctype = char (ones (1,n*J)*('B'));
   
   cplex.addCols(obj, [], lb, ub, ctype);
   
   % Each Disk must be assigned to exactly one speed
   for i = 1:n
      speedLimit = zeros (1,J*n);
      speedLimit(((i-1)*J+1):(i*J)) = ones (1, J);
      cplex.addRows(1, speedLimit, 1);
   end
   
   % The average response time should not be more than Rlimit
   for i = 1:J
      v = zeros (1, J*n);
      v(i:J:end) = TotalReq(:,i);
      cplex.addRows(-inf, v, TotalN*Rlimit);
   end
   
   cplex.solve();
   cplex.writeModel('Hibernator.lp');
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   if cplex.Solution.status == 101
      fprintf ('\nMinimum Energy: %f \n', cplex.Solution.objval);
      EnergySolution = reshape (cplex.Solution.x, J, n);
   end
catch m
   throw (m);
end
