function res = bc_bvp_cat(ya,yb,beta)
global xp
ds = 1e-4;
res=[ya(5) ;
    ya(1)-xp;
    yb(5);
    ya(2);
    ya(3)-pi/2;
    yb(1)-ds
    yb(3)-pi
    ];
% res=[yb(3) - pi;
%     ya(1)-ds;
%     yb(2) - 1;
%     ya(2);
%     ya(3)-0;
%     yb(1)-ds
%     yb(6)-0
%     ];
end
