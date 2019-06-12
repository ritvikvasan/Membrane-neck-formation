function rhs = rhs_bvp(x,y, beta)
global c0 a0 g delta p pw wd r0 z0 f0 aF k0 aF_neck

% spontaneous curvature
%c = -0.5*c0*(tanh(g*(x+a0))-tanh(g*(x-a0)));
%c = -0.5*c0*(1 - tanh(g*(flip(x) - a0))); 
%c = -0.5*c0*(1 - tanh(g*(x - a0)));
%c = -0.5*c0*(tanh(g*(x - z0/r0*0.5 + a0)))+0.5*c0*(tanh(g*(x - z0/r0*0.5- a0))) ; 
c = -c0/2 * ( 1./(1+exp((x-a0)/0.005)));

% derivative of spontaneous curvature
%dc = -0.5*c0*g*(1 - tanh(g*(x+a0))^2) + 0.5*c0*g*(1-tanh(g*(x-a0))^2);
%dc = 0.5*c0*g*(tanh(g*(flip(x) - a0))^2 - 1);
%dc = -0.5*c0*g*(tanh(g*(x - a0))^2 - 1);
%dc = -0.5*c0*g*(1-tanh(g*(x - z0/r0*0.5 + a0))^2)+ 0.5*c0*g*(1 -tanh(g*(x - z0/r0*0.5 - a0))^2);
dc = -c0/2 * (-200*(exp(200*a0 + 200*x))./(exp(200*a0)+exp(200*x)).^2);

% bending modulus
b = 1 + 0.5*(delta-1)*(1 - tanh(g*(x - a0)));
%b = k;
% opposing pressure from rigid wall
dp = 0; % -p*(tanh(g*(y(2)))) - pw/2*(1 + tanh(g*(y(2) - wd)));

%force 
%fbar = -0.5*beta/aF*(tanh(g*(x - 10/5*0.5 + aF)))+0.5*beta/aF*(tanh(g*(x - 10/5*0.5- aF)));
%fbar = -beta*(0.5*(1 - tanh(g*(x - aF))))/aF;
%fbar = -beta*(0.5*(1 - tanh(g*(x - aF))))/aF;
% fbar = -beta/(2*aF) * ( 1./(1+exp((x-aF)/0.005)));
fbar = -beta/aF_neck * ( 1./(1+exp((x-aF_neck)/0.005)));

rhs=[cos(y(3))/y(1);
     sin(y(3))/y(1);
     2*y(4)/y(1) - sin(y(3))/y(1)^2;
     y(5)/y(1)^2 + dc;
     dp/b - sin(y(3))*fbar + 2*y(4)*((y(4)-c)^2+y(6)/b)-2*(y(4)-c)*(y(4)^2+(y(4)-sin(y(3))/y(1))^2);
     2*b*(y(4)-c)*dc - cos(y(3))*fbar/y(1)
     ];
end