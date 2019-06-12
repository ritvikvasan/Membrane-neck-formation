function in = init_bvp_cat(x)
global z0 r0
m = length(x);
in = [1/cos(0)*cos((x)).*ones(1,m);
    sin((x)).*ones(1,m);
    pi/2*ones(1,m);
    ones(1,m);
    zeros(1,m);
    1/4*ones(1,m)
    ];
end