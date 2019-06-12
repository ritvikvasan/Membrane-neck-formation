function res = bc_bvp_2(ya,yb)
global z0 r0
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

% this is correct
res=[ya(5);
    yb(3)-pi/2;
    yb(1)-1;
    ya(2);
    ya(3)-pi/2;
    yb(2)-z0/r0
    ];

% for haleh/ full sphere
% res=[ya(1) - ds;
%     yb(3)-pi;
%     yb(1)-ds;
%     ya(2);
%     ya(3)-0;
%     yb(6) - 0
%     ];
end
