bp = 0.95;
nums = 55;
z_c = zeros(1,length(XP)); 
for ii = 1:nums
    temp = loopsol(8,:,ii)*20;    
    meshes = linspace(bp, 2, 1001);
    z_c1 = temp(find(meshes > bp + 0.05, 1));
    z_c2 = min(temp);
    z_c(ii) = z_c1 - z_c2;
end