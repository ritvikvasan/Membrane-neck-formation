function [W, E, Epos,Epos2,Epos3, H, r, lam] = plotenergy(Sol , pos, alpha, l, XP)

R0=20;
k0 = 320;
%k2 = -0.53*k0;
k2=0;
%mesh = (0:1/2000:1).^2;
mesh = 0:0.001:1;
t = alpha*mesh;
H(:,:) = Sol(4,:,:)/R0;
r(:,:) = Sol(1,:,:)*R0;
y(:,:) = Sol(2,:,:)*R0;
lam(:,:) = Sol(6,:,:)*k0/R0^2;
s(:,:) = Sol(3,:,:);
tt = t;

for j = 1:pos
    for i = 1:1000
    K(i,j) = H(i, j)^2 - (H(i,j) - sin(s(i,j))/r(i,j))^2;
    W(i,j) = k0*(H(i,j))^2 + k2*K(i,j) ; 
    %     E(i,j) = W(i,j)*2*pi*r(i,j)*abs((y(i+1,j) - y(i,j))) ;
    E(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2 + lam(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    E2(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    E3(i,j) = lam(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
    end
    Epos(j) = sum(E(:,j));
    Epos2(j) = sum(E2(:,j));
    Epos3(j) = sum(E3(:,j));
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
% plot(XP/20, Epos*20/(2*pi*20*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('x')
% ylabel('E/(2\pik) (including tension)')
% axis square
% hold off
% 
% 
% fighandle = figure(2);
% hold on
% plot(XP/20, Epos3*20/(2*pi*20*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('x')
% ylabel('E/(2\pik) (only tension)')
% axis square
% hold off
%xlim([20 70])
% 
% 
% fighandle = figure(2);
% hold on
% plot(XP/20, l*20/320)
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('x')
% ylabel('FR/k')
% axis square
% hold off
% 
% 
% 
% 
% 
% 
% 
% 
% fighandle = figure(5);
% hold on
% plot(XP/20, Epos2*20/(2*pi*20*320)) % E = 2977.2 kT (for breakt at 7, tube Z=800 ithout force
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',14, 'fontweight','bold')
% xlabel('x')
% ylabel('E/(2\pik)')
% axis square
% hold off
% 
% 



end

