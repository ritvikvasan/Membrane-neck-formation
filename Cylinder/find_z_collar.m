z_c = zeros(1,length(XP)); 
for ii = 1:123
    temp = loopsol(2,:,ii)*20;    
    meshes = (0:0.001:1)*2;
    z_c(ii) = temp(find(meshes == 0.05, 1));
end