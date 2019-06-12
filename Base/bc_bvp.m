function res = bc_bvp(ya,yb, beta)
global xp z0 r0 lam
ds = 1e-4;
%{
res=[yb(2)-z0/r0;
    ya(3)-pi/2;
    ya(2);
    yb(1) - xp;
    yb(5);
    yb(3)-pi/2;
    ya(6)-1/4
    ];
%}
% res=[ya(5);
%     yb(9)-pi/2;
%     ya(2);
%     ya(7) - xp;
%     yb(8)-z0/r0;
%     ya(3)-pi/2;
%     yb(12)-0.25;%yb(7)-1;
%     yb(1) - ya(7);
%     yb(2) - ya(8);
%     yb(3) - ya(9);
%     yb(4) - ya(10);
%     yb(5) - ya(11);
%     yb(6) - ya(12);
%     ];

% this is correct
% res=[ya(5);
%     yb(9)-0;
%     ya(2);
%     ya(7) - xp;
%     yb(7)-1.5;%yb(8)-z0/r0;
%     ya(3)-pi/2;
%     yb(12)-lam;%yb(7)-1;
%     yb(1) - ya(7);
%     yb(2) - ya(8);
%     yb(3) - ya(9);
%     yb(4) - ya(10);
%     yb(5) - ya(11);
%     yb(6) - ya(12);
%     ];
% this is test im an idiot this is really correct
 res=[yb(11);
    ya(3)-pi;
    ya(2);
    yb(1) - xp;
    ya(1)-1.5;%yb(8)-z0/r0;
    yb(9)-pi/2;
    ya(6)-lam;%yb(7)-1;
    yb(1) - ya(7);
    yb(2) - ya(8);
    yb(3) - ya(9);
    yb(4) - ya(10);
    yb(5) - ya(11);
    yb(6) - ya(12);
    ];
end
