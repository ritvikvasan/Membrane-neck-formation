function in = init_bvp(x)
global aa
m = length(x);
tt = (2 - aa)/aa;

in = [ones(1,m);
    x;
    pi/2*ones(1,m);
    1/2*ones(1,m);
    zeros(1,m);
    1/4*ones(1,m);
    ones(1,m);
    tt*(x) + aa;
    pi/2*ones(1,m);
    1/2*ones(1,m);
    zeros(1,m);
    1/4*ones(1,m)
    ];
end