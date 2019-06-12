function res = bc_bvp(ya,yb, beta)
global xp z0 r0
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
res=[ya(5);
    yb(3)-pi/2;
    ya(2);
    ya(1) - xp;
    yb(2)-z0/r0;
    ya(3)-pi/2;
    yb(1)-1
    ];
end
