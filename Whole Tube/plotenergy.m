function [W, E, Epos,Epos2,Epos3,Epos4, H, r, lam,I,LT, testf, testf2, Elamt] = plotenergy(Sol , pos, alpha, FvsZp)

R0=4;
k0 = 320;
%k2 = -0.53*k0;
k2=0;
mesh = (0:1/2000:1).^2;
t = alpha*mesh;

H(:,:) = [Sol(4,:,:) Sol(10,:,:)]/R0;
r(:,:) = [Sol(1,:,:) Sol(7,:,:)]*R0;
y(:,:) = [Sol(2,:,:) Sol(8,:,:)]*R0;
lam(:,:) = [Sol(6,:,:) Sol(12,:,:)]*k0/R0^2;
s(:,:) = [Sol(3,:,:) Sol(9,:,:)];
c = 0;
tt = [t (200-alpha)/alpha*t + alpha];
syms aa


for j = 1:pos
    E4(j) = FvsZp(2,j)*800;
    for i = 1:4001
    K(i,j) = H(i, j)^2 - (H(i,j) - sin(s(i,j))/r(i,j))^2;
    W(i,j) = k0*(H(i,j))^2 + k2*K(i,j) ; 
    %     E(i,j) = W(i,j)*2*pi*r(i,j)*abs((y(i+1,j) - y(i,j))) ;
    E(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2 + lam(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    E2(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    E3(i,j) = lam(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    end
    psip1(j) = 2*r(2000,j)*H(2000,j) - sin(psi(2000,j))./r(2000,j);
    psip2(j) = 2*r(2002,j)*H(2002,j) - sin(psi(2002,j))./r(2002,j);
%     psip1(j) = sin(psi(2000,j))./r(2000,j);
%     psip2(j) = sin(psi(2002,j))./r(2002,j);
    %testf(j) = (psip1(j) + psip2(j) - 1)./(lam(2002,j) + 1/2*(2*psip1(j) - 1)^2);
    %testf(j) = 2/20*(lam(2001,j) - 1/(r(2001,j))^2);
    %testf(j) = 2/20 - 1/(r(2000,j)) - 1/(r(2002,j));
    %testf(j) = 2*mean(lam(1965:2001,j))*R0^2/(k0*(70/R0 + r(2001,j)/R0)) - 2*(70/R0)^2/((r(2001,j)/R0)^2*(70/R0+r(2001,j)/R0)^3);
    %testf(j) = 2*mean(lam(1965:2001,j))/((70 + r(2001,j))) - 2*250^2/((r(2001,j)/R0)^2*(250/R0+r(2001,j)/R0)^3);
    testf(j) = 2*lam(1980,j)/(r(2001,j) + 1) - 2*lam(1980,j)/(r(2001,j)^2*(1+r(2001,j))^3); 
    Epos(j) = sum(E(:,j)) - E4(j);
    Epos2(j) = sum(E2(:,j));
    I(j) = find(t-alpha+5>=0,1);
    Epos4(j) = sum(E(I(j):2001,j));
    Epos3(j) = sum(E3(:,j));
    LT(j) = (H(2001,j) + H(2002,j))*(2*320);
    testf2(j) = -testf(j)/wrightOmega(-log(-1/testf(j)));
    for ii = 1965:2001
        Elam(ii,j) = lam(ii,j)*(tt(ii+1) - tt(ii))*2*pi*R0^2;
    end
    Elamt(j) = sum(Elam(:,j));
%     Epos2(j) = sum(E(1901:21011,j ));
end

% for j = 1:pos
%     for i = 1:4001
%     K(i,j) = H(i, j)^2 - (H(i,j) - sin(s(i,j))/r(i,j))^2;
%     W(i,j) = k0*H(i,j)^2 + k2*K(i,j);
%     E(i,j) = W(i,j)*2*pi*r(i,j)*abs((y(i+1,j) - y(i,j)));
%     end
%     Epos(j) = sum(E(:,j));
%     Epos2(j) = sum(E(1901:21011,j ));
% end


% fighandle = figure(1);
% hold on
% plot(FvsZp(1,:)/20, Epos*20/(2*pi*800*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('E/E(1) (total energy)')
% axis square
% hold off
%xlim([20 70])

% 
% fighandle = figure(2);
% hold on
% plot(FvsZp(1,:)/20, FvsZp(2,:)*20/320)
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('Axial force (pN)')
% axis square
% hold off
% %xlim([20 70])
% 
% 
% fighandle = figure(3);
% hold on
% plot(FvsZp(1,:)/20, FvsZp(3,:)*20/320)
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('Radial force (pN)')
% axis square
% hold off
% %xlim([20 70])
% 
% 
% % fighandle = figure(4);
% % hold on
% % plot(FvsZp(1,:), sqrt(FvsZp(2,:).^2+ FvsZp(3,:).^2),'.')
% % set(fighandle, 'Position', [0, 1000, 300, 300]);
% % set(gca, 'fontsize',14, 'fontweight','bold')
% % xlabel('Radius (nm)')
% % ylabel('Total force (pN)')
% % axis square
% % hold off
% 
% %xlim([20 70])
% 
% 
% fighandle = figure(15);
% hold on
% plot(FvsZp(1,:)/20, Epos3*20/(2*pi*800*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('Only tension')
% axis square
% hold off
%xlim([20 70])

% fighandle = figure(5);
% hold on
% plot(FvsZp(1,:)/20, Epos2*20/(2*pi*800*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('E2/E2(1) (only bending)')
% axis square
% hold off
%xlim([20 70])



% 
% fighandle = figure(6);
% hold on
% plot(FvsZp(1,:)/20, Epos4*20/(2*pi*800*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('E4/E4(1) (only b energy)')
% axis square
% hold off
% %xlim([20 70])
% 
% fighandle = figure(7);
% hold on
% plot(FvsZp(1,:)/20, LT) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('Radius (nm)')
% ylabel('LT (average mean curv)')
% axis square
% hold off


end

