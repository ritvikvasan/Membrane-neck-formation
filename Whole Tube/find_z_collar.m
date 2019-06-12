bp = 130;
num = 63;
z_c = zeros(1,num); 
for ii = 1:num
    temp = coatPullSol(2,:,ii)*4;    
    meshes = bp*(0:1/2000:1).^2;
    z_c1 = temp(find(meshes > 1, 1));
    z_c2 = temp(1);
    z_c(ii) = z_c2 - z_c1;
end