function in = init_bvp(x)
m = length(x);
in = [ones(1,m);
    x;
    pi/2*ones(1,m);
    1/2*ones(1,m);
    zeros(1,m);
    1/4*ones(1,m)
    ];
end