%% File to run homework 1 deliverables

clc; clear; 

format longg

mu_sun = 1.327124400189*10^11; % in km^3/s^2
au = 1.49597870691*10^8; % km
DU = au;

P = 2*pi*sqrt(au^3/mu_sun); % period
TU = sqrt(au^3/mu_sun); % seconds

mu_sun = 1; % use with DU and TU computed above

% Units are au and au/day
% Initial state vectors for epoch 2017-Jan-01 00:00:00.0000 UTC
r1I = [3.515868886595499 * 10^-2, -3.162046390773074, 4.493983111703389];
v1I = [-2.317577766980901 * 10^-3, 9.843360903693031 * 10^-3, -1.541856855538041 * 10^-2];
r2I = [7.249472033259724, 14.61063037906177, 14.24274452216359];
v2I = [-8.241709369476881 * 10^-3, -1.156219024581502 * 10^-2, -1.317135977481448 * 10^-2];

% Initial state vector for the Earth
rE = [-1.796136509111975 * 10^-1, 9.667949206859814 * 10^-1, -3.668681017942158 * 10^-5];
vE = [-1.720038360888334 * 10^-2, -3.211186197806460 * 10^-3, 7.927736735960840 * 10^-7];

% Convert to canonical units (only need to convert velocity)
v1I = (v1I/86400)*TU;
v2I = (v2I/86400)*TU;
vE = (vE/86400)*TU;

%% Homework 1 Problem 3
% Propagate Earth departure from January 2017 - December 2017
% Propagate 1I arrival from August 2017 - January 2019
% Solve Lambert's problem from all departure to arrival dates
% Calculate DV from Lambert's problem
% Each departure arrival pair will have an associated DV
% Make pork chop plot
% 
% % Testing Earth propagation around sun for 1I

dt_e1 = datetime(2017,1,1):days(1):datetime(2017,12,31);

r_traj_e1 = zeros(length(dt_e1), 3);
v_traj_e1 = zeros(length(dt_e1), 3);

for i=1:length(dt_e1)
    tof = days(dt_e1(i)-dt_e1(1))/(TU/86400);
    [r_traj_e1(i,:),v_traj_e1(i,:)] = r_v(mu_sun, tof, rE, vE);
end

% Testing 1I propagation wrt to sun

dt_1I = datetime(2017,8,1):days(1):datetime(2019,1,1);

r_traj_1I = zeros(length(dt_1I), 3);
v_traj_1I = zeros(length(dt_1I), 3);

for i=1:length(dt_1I)
    tof = days(dt_1I(i)-dt_e1(1))/(TU/86400);
    [r_traj_1I(i,:),v_traj_1I(i,:)] = r_v(mu_sun, tof, r1I, v1I);
end

dv_1I = zeros(length(dt_1I), length(dt_e1));

for i=1:length(dt_1I) % arrival
    for j=1:length(dt_e1) % departure
        [i,j]
        TOF = days(dt_1I(i)-dt_e1(j))/(TU/86400);
        if TOF > 0
            traj_type = 'prograde';
            [v1,v2] = lambert_solver(mu_sun, r_traj_e1(j,:), r_traj_1I(i,:), TOF, traj_type);
            if ~isnan(v1(1)) && ~isnan(v2(1))
                dv_1I(i,j) = norm(v1 - v_traj_e1(j,:));% + norm(v_traj_1I(i,:) - v2);
            end
        end
    end
end

%% Problem 4

% % Earth propagation for 2I
% 
% dt_e2 = datetime(2017,1,1):days(1):datetime(2020,7,1);
% 
% r_traj_e2 = zeros(length(dt_e2), 3);
% v_traj_e2 = zeros(length(dt_e2), 3);
% 
% for i=1:length(dt_e2)
%     tof = days(dt_e2(i)-dt_e2(1))/(TU/86400);
%     [r_traj_e2(i,:),v_traj_e2(i,:)] = r_v(mu_sun, tof, rE, vE);
% end
% 
% % Testing 2I propagation wrt to sun
% 
% dt_2I = datetime(2019,6,1):days(1):datetime(2022,1,1);
% 
% r_traj_2I = zeros(length(dt_2I), 3);
% v_traj_2I = zeros(length(dt_2I), 3);
% 
% for i=1:length(dt_2I)
%     tof = days(dt_2I(i)-dt_e2(1))/(TU/86400);
%     [r_traj_2I(i,:),v_traj_2I(i,:)] = r_v(mu_sun, tof, r2I, v2I);
% end
% 
% dv_2I = zeros(length(dt_2I), length(dt_e2));
% 
% for i=1:length(dt_2I) % arrival
%     for j=1:length(dt_e2) % departure
%         [i,j]
%         TOF = days(dt_2I(i)-dt_e2(j))/(TU/86400);
%         if TOF > 0
%             traj_type = 'prograde';
%             [v1,v2] = lambert_solver(mu_sun, r_traj_e2(j,:), r_traj_2I(i,:), TOF, traj_type);
%             if ~isnan(v1(1)) && ~isnan(v2(1))
%                 dv_2I(i,j) = norm(v1 - v_traj_e2(j,:));% + norm(v_traj_2I(i,:) - v2);
%             end
%         end
%     end
% end

% Plotting test sanity check
% figure(1)
% hold on 
% grid on
% axis equal
% % plot3(r_traj_e1(:,1), r_traj_e1(:,2), r_traj_e1(:,3))
% % plot3(r_traj_e1(end,1), r_traj_e1(end,2), r_traj_e1(end,3),'.','MarkerSize',10)
% plot3(0,0,0,'.','MarkerSize',20)
% % plot3(r_traj_1I(:,1), r_traj_1I(:,2), r_traj_1I(:,3))
% % plot3(r_traj_1I(end,1), r_traj_1I(end,2), r_traj_1I(end,3), '.','MarkerSize',10)
% % % plot3(r_traj_e1(j,1), r_traj_e1(j,2), r_traj_e1(j,3), '.', 'MarkerSize',5)
% % % plot3(r_traj_1I(i,1), r_traj_1I(i,2), r_traj_1I(i,3), '.', 'MarkerSize',5)
% plot3(r_traj_e2(:,1), r_traj_e2(:,2), r_traj_e2(:,3))
% plot3(r_traj_2I(:,1), r_traj_2I(:,2), r_traj_2I(:,3))
% % plot3(0,0,0,'.','MarkerSize',10)
% % 

% figure(1)
% dv_kms = dv*DU/TU;
% % dv_kms = dv_kms(dv_kms<50);
% dv_kms(dv_kms>20) = 0;
% contourf(dv_kms,'edgecolor','none');
% colorbar
% % set(gca,'XLim',[1 365])
% % set(gca, 'YLim', [1,760-212])
% 
figure(4)
dv_kms = dv_1I*DU/TU;
dv_kms(dv_kms>20) = NaN;
idx = find(dv_kms==0);
dv_kms(idx) = NaN;
surf(dt_e1, dt_1I, dv_kms, 'EdgeColor', 'Interp', 'FaceColor', 'Interp')
colormap('jet')
colorbar
view([0 90])
% set(gca,'XLim',[1 365])
% set(gca, 'YLim', [1,760-212])

%xtickformat("yyyy-MM-dd")

% s_date = datenum('01-01-2017');
% e_date = datenum('12-31-2017');

% dateaxis('x', 12, datetime(2017,1,1))

% datetick('x',3)

% datetick('y',3)

% dateaxis('x',3,datetime(2017,1,1))
% xticks(months(['jan', 'feb']))

% figure(2)
% contour(dv)
% colorbar
