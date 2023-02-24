clc; clear; 

load('hw_1.mat', 'dv')

% [r,c,~] = size(v1_lambert);
% dv_1 = zeros(r, c, 1);
% dv_2 = zeros(r, c, 1);
% for i=1:r
%     for j=1:c
%         tmp = v1_lambert(i,j,:);
%         if norm([tmp(1) tmp(2) tmp(3)]) == 0
%             dv1_tmp = v1_lambert(i,j,:) - v_traj_e1(i);
%             dv_1(i,j,:) = norm([dv1_tmp(1) dv1_tmp(2) dv1_tmp(3)]);
%             dv2_tmp = v_traj_1I(i) - v2_lambert(i,j,:);
%             dv_2(i,j,:) = norm([dv2_tmp(1) dv2_tmp(2) dv2_tmp(3)]);
%         end
%     end
% end

% dv = v1_lambert + v2_lambert;

val = min(min(dv));

mu_sun = 1.327124400189*10^11; % in km^3/s^2
au = 1.495978707*10^8; % km
DU = au;

P = 2*pi*sqrt(au^3/mu_sun); % period
TU = sqrt(au^3/mu_sun); % seconds

% dv = dv*DU/TU;
% dv(dv>50) = 0;
% dv = dv(dv<50);

idx = find(dv>100)
dv(idx) = 0;

figure(1)
colormap(hot)
contourf(dv)
colorbar
dateaxis('x',1,datetime(2017,1,1))
dateaxis('y', 1, datetime(2017, 8,1))

