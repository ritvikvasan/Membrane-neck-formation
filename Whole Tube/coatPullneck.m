% Inputs: 
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   alpha0 - dimensionless coat area
%   rF - radius of the area of applied force (on a flat membrane), in units of nm
%   zp - pole z-position, in units of nm
%   f0 - Initial guess for the applied force, in units of pN
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the solution

% Outputs:
%   t - area mesh points
%   Sol - solution array
%   f - force to hold pole at specified z-position, in units of pN

function [t,Sol,faxial, fradial, Z_c] = coatPullneck(alpha, mesh, lambda, alpha0, rF, XP, f0, k0, dk, P, gamma, C0, R0, initSol, Fz, zp)

t=alpha*mesh;   % area mesh points

% declare and assign global variables to be used in nested functions
global a gt a0 iSol lam c0 g delta p aF xp fz tt yp aTop aF_neck

a = alpha;
a0 = alpha0;
gt = t;
lam = lambda*R0^2/k0;
c0 = R0*C0;
g = gamma;
delta = dk;
p = P*R0^3/k0;
aF = ((rF/R0)^2)/2;
xp = XP/R0;
iSol = initSol;
fz = Fz*R0^3/(k0*2*pi*R0^2);
yp = zp/R0;

% Find Z_collar
temp = initSol(2,:)*R0;
Z_c = abs(temp(end) - temp(find(t>=alpha-1, 1)));
% disp('Height of collar is', 1)
aF_neck = a - t(find(temp<=min(temp)+1, 1));
aF_neck = 1;
disp(aF_neck)
% Z_c = 1;
disp(Z_c)
% beta = [f0*R0^3/k0 f0*R0^2*Z_c/k0];
beta = [f0*R0^3/k0 f0*R0^3/k0];
% beta = f0*(R0^3/k0)/((Afd)*2*pi*R0^2);
tt = (200 - alpha)/alpha;
%aTop = t(find(initSol(2,:)*R0>=(initSol(2,end)*R0+10),1,'last'));
%tt = (800 - alpha)/alpha;

% initial guess structure
solinit = bvpinit(t,@mat4init, beta);

% solver options; increasing maximum number of mesh points solver will use
options = bvpset('NMax', 100*length(t), 'RelTol', 1e-3);

% solve the boundary value problem
sol = bvp4c(@mat4ode,@mat4bc,solinit,options);

% extract force
faxial = sol.parameters(1)*k0/R0^3*2*pi*R0^2;
fradial = sol.parameters(2)*k0/R0^3;
%f = sol.parameters*k0/R0^3*pi*rF^2;

% evaluate the solution on the original mesh points
Sol = deval(sol,t);

% plot the resultant profile of the membrane
coatArea = [200/(tt+1)-5 200/(tt + 1)];
%coatArea = [aTop 200/(tt + 1)];
%coatArea = [800/(tt+1)-20 800/(tt + 1)];
actArea = [0 aF];
plotTitle = sprintf('Membrane profile, f = %0.3f pN', fradial);
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
xLim = [-400 400];
% yLim = [0 600];

plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], xLim, [], plotTitle, 0)


%%Define the variables
% X(1)=x;
% X(2)=y;
% X(3)=phi;
% X(4)=H;
% X(5)=L;

%%%%%the differential equations
%------------------------------------
function dXdt = mat4ode(t, X, beta)
%parameters
global c0 a0 g delta p aF fz tt aTop aF_neck

% spontaneous curvature
%c = 0.5*c0*(1 - tanh(g*(t - a0))); 
%c = -beta(1);
c = -c0;

% derivative of spontaneous curvature
%dc = 0.5*c0*g*(tanh(g*(t - a0))^2 - 1);
dc = 0;

% bending modulus
b = 1 + 0.5*(delta-1)*(1 - tanh(g*(t - a0)));

% opposing pressure from rigid wall
% this is correct
dp = -p*(tanh(g*(X(2))));

% test
%dp = -beta(2)*(tanh(g*(X(2))));
%dp = -p*(tanh(g*(X(2)))) - pw/2*(1 + tanh(g*(X(2) - wd)));

% no wall to oppose pressure
%dp = p; 

% applied force
fbar = beta(1)*(0.5*(1 - tanh(g*(t - aF))))/aF;
%fbar = fz*(0.5*(1 - tanh(g*(t - aF))))/aF;

aF2 = 200./(tt + 1);
%aF2 = 800./(tt + 1);
%aBet = aF2 - aTop;

% fbar2 = -beta*(0.5*(1 - tanh(g*(t - 40 - aF))))/(40+aF);
%fbar2 = -beta(2)*(0.5*tanh(g*(t - (aF2-aBet))) - 0.5*tanh(g*(t - aF2 - 0.1)));

% this is correct
%fbar2 = -beta(2)*(0.5*tanh(g*(t - (aF2-5))) - 0.5*tanh(g*(t - aF2 - 0.1)));
%fbar2 = -beta(2)*(0.5*tanh(g*(t - (aF2-1))) - 0.5*tanh(g*(t - aF2 - 0.1)));
%fbar = -beta * ( 1./(1+exp((x-aF)/0.005)));
fbar2 = beta(2)/aF_neck *( ( 1./(1+exp((t-(aF2-aF_neck))/0.005))) - ( 1./(1+exp((t-aF2)/0.005))));
% test
%fbar2 = 0;
% fbar3 = -beta*(0.5*(1 - tanh(g*(t - 0.5))));
fbar3 = 0;

dXdt = [cos(X(3))/X(1)
        sin(X(3))/X(1)
        (2*X(1)*X(4) - sin(X(3)))/X(1)^2
        X(5)/X(1)^2 + dc
        dp/b + fbar*cos(X(3)) - fbar2*sin(X(3)) + 2*X(4)*((X(4)-c)^2 + X(6)/b) - 2*(X(4)-c)*(X(4)^2 + (X(4)-sin(X(3))/X(1))^2)
        2*b*(X(4) - c)*dc - fbar*sin(X(3))/X(1) - fbar2*cos(X(3))/X(1)
        tt*cos(X(9))/X(7)  %3 times with 3t + 50, 4 times with 4t + 40
        tt*sin(X(9))/X(7)
        tt*(2*X(7)*X(10) - sin(X(9)))/X(7)^2
        tt*(X(11)/X(7)^2 + dc)
        tt*(dp/b -fbar3*sin(X(9)) + 2*X(10)*((X(10))^2 + X(12)/b) - 2*(X(10))*(X(10)^2 + (X(10)-sin(X(9))/X(7))^2))
        tt*(2*b*(X(10))*dc - fbar3*cos(X(9))/X(7) )
         ];
            

%-------------------------boundary conditions-------------

function res = mat4bc(Xa,Xb,beta) % y position at pole specified
global lam xp yp
ds = 1e-4;
   
    res = [ Xa(1) - ds
            Xb(1) - xp%Xa(2) - 30
            Xb(8)
            Xb(9)
            Xa(3)               
            Xa(2) - yp %Xb(9)             
            Xa(5) 
            Xb(12) - lam
            Xb(1) - Xa(7)
            Xb(2) - Xa(8)
            Xb(3) - Xa(9)
            Xb(4) - Xa(10)
            Xb(5) - Xa(11)
            Xb(6) - Xa(12)
            ];
        
        
%-----------------------------------Initial guesses------------


function Xinit = mat4init(t)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Xinit = iSol(:,find(gt==t));