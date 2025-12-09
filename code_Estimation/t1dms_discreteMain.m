close all
clear;
clc;

Path = [
    './t1dms_data/subjects/:',...
    './t1dms_data/simu_data/:'];

% Path = [
%     '.\t1dms_data\subjects\;',...
%     '.\t1dms_data\simu_data\;'];

addpath(Path);

table = [];
%%%%
feedback = "uva";
for subject_num =3 %subject
    for seed = 9 %seed

        d = "ssogmm"; %ssogmm inequality linear


        if subject_num == 10
            load(['adult#0',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        else
            load(['adult#00',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        end

        load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])

        %seed9によるパラメータ推定結果を利用
         % load(['t1dms_adult',num2str(subject_num),'_seed9.mat'],'p_esti','x_esti')


        %%%
        Ts = 5;
        n_day = 4;
        t_end = 1440*n_day;

        %data procssing
        % modes = data_ssogmm.mode;
        if Ts==1
            modes = repelem(modes,5);
        elseif Ts==2.5
            modes = repelem(modes,2);
        end

        N=size(modes,2);
        t_span = 200:Ts:200+1075*5-1;
        N_step = 1:N;

        %%%
        G = data.G.signals.values';

        insulin = data.injection.signals.values';
        meal = 1000*data.CHO.signals.values;

        %%%   ssogmm
        model = ssogmmCL("u_inb",(1/6)*Ib,...
            "BW",BW,...
            "G_b",Gb);

        model.p_0 = p_esti;
        model.x_0 = [x_esti(:,1)];%

        data_ssogmm = simu_t1dm_ssogmm(Ts,t_end,t_span,modes,...
            "insulin",insulin,...
            "meal",meal,...
            "ssogmm",model,"dynamics",d,"q",q_esti);
        ssogmm = data_ssogmm.xs([1 2 5 6 7],:);

        %%%%%

        %discrete simulation
        %%%common setting
        Xs = data_ssogmm.xs;

        %estimate Xinit
        x_ini = x_esti(:,500);
        x_ini = x_ini([1 2 5 6 7]);

        %feedback
        if feedback == "uva"
        Y = G(data_ssogmm.ts);
        elseif feedback == "ssogmm"
        Y = data_ssogmm.Gs;
        end

        %%input
        p = p_esti;
        u_m = data_ssogmm.meal;
        u_i = data_ssogmm.insulin;
        modes = data_ssogmm.modes;

        %%Ra
        R_a = get_Ra(Xs(4,:),p,modes,N);

        [xe,A_d,B_d1,B_d2,C,i_b,Ge] = linear_matrix(p,Ts);
        U = [u_m;u_i;R_a;modes];
        

        %%%

        %%%linear
        delta_x=zeros(5,N);
        x_linear=zeros(5,N);

        delta_x(:,1) = Xs([1 2 5 6 7],1)-xe;
        x_linear(:,1) = Xs([1 2 5 6 7],1);

        for k=2:N

            delta_u_i = u_i(:,k-1)-i_b;

            delta_x(:,k) = A_d*delta_x(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1);
            x_linear(:,k) = delta_x(:,k) + xe;
        end


        %observer
        % 状態推定誤差の重み
        Q = diag([1, 1e-6, 1e5, 1e5, 1e3]);%
        R = 1;
        [L,pole,pole_ori] = gain(p,Ts,Q,R);

        delta_xob = zeros(5,N);
        xob = zeros(5,N);

        delta_xob(:,1) = x_ini - xe;
        xob(:,1) = x_ini;

        delta_y = Y - Ge;
        for k=2:N

            delta_u_i = u_i(:,k-1) - i_b;

            delta_xob(:,k) = A_d*delta_xob(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1) + L*(delta_y(:,k-1) - C*delta_xob(:,k-1));
            xob(:,k) = delta_xob(:,k) + xe;
        end



        %%ssogmm
        dic = dic_ssogmmCL(p,Ts,N,U);
        x_ssogmm = dic.dicrete(Xs,"ZOH");


        %kalman
        A = dic.A;
        A_di = dic.A_di;
        A_dx = dic.A_dx;
        C = [1 0 0 0 0];

        xhat = zeros(5,N);

        a = [10,1e-6,1000,1000,500];
        P = diag(a);
        xhat(:,1) = x_ini;

        %%
        Q = 1e-2; %システム雑音
        R = 0.006; %観測雑音
        % v = randn(1,N)*sqrtm(Q);
        % w = randn(1,N)*sqrtm(R);

        B = [0.1, 1e-6, 100, 100, 25];
        B_diag = diag(B);
        g_block = [];
        v = [];
        v = [v,Y(1,1)-x_ini(1,1)];


        NIS = zeros(1,N);
        NEES= zeros(1,N);
        for k = 2:N

            xhatm = dic.dynamics(xhat(:,k-1),k);
            v = [v,Y(:,k) - xhatm(1,1)];

            A_tmp = [zeros(3,2) A_di];
            A_temp =[A(xhat(1,k-1),xhat(2,k-1)) ; A_tmp];
            % Pm = A(xhat(1,k-1),xhat(2,k-1))*P*A(xhat(1,k-1),xhat(2,k-1)) + b*Q*b';
            Pm = A_temp*P*A_temp' + B_diag*Q*B_diag';

            nu = Y(:,k)-C*xhatm;
            S = C*Pm*C' + R;
            NIS(k) = nu'*(S\nu);

            g = Pm*C'/(C*Pm*C' + R);

            xhat(:,k) = xhatm + g*(Y(:,k) - C*xhatm);
            P = (eye(5) - g*C)*Pm;

            e = ssogmm(:,k) - xhat(:,k);
            NEES(k) = e'*(P\e);
            g_block = [g_block,g];
        end
 
    end
    obs = rmse(xob(1,:)',Y(1,:)');
ekf = rmse(xhat(1,:)',Y(1,:)');

d = floor(N - N/4);
obs4 = rmse(xob(1,d:end)',Y(1,d:end)');
ekf4 = rmse(xhat(1,d:end)',Y(1,d:end)');
table = [table,[obs;obs4;ekf;ekf4]];
end


% mean(NIS)
% mean(NEES)

for i=1:5
    figure
    plot(data_ssogmm.ts,ssogmm(i,:),'LineWidth',4)
    hold on
    % plot(data_ssogmm.ts,xob(i,:),'.')
    plot(data_ssogmm.ts,xhat(i,:),'.','MarkerSize',13)
    h_axes = gca;
h_axes.XAxis.FontSize = 25;
h_axes.YAxis.FontSize = 25;

if i == 1
ylim([110 190]);
elseif i == 2
ylim([-1*10^-3 7*10^-3]);
elseif i == 3
ylim([3000 9000]);
elseif i == 4
ylim([3000 6500]);
elseif i == 5
ylim([450 1000]);
end

end




% figure
% tiledlayout(2,3)
% for i=1:5
%     nexttile
% stairs(g_block(i,:))
% end

for i = 1:5
    figure(i)
    if i == 1
        hold on
        plot(G,'LineWidth',4)
        legend({'SSOGMM','EKF','UVA/Padova'},'FontSize',15,'FontWeight','bold')
    else
        legend({'SSOGMM','EKF'},'FontSize',15,'FontWeight','bold')
    end
    grid on
    xlim([0 3400])
    xlabel('Time[min]')
    if i == 1
        ylabel('Blood Glucose[mg/dl]')
    elseif i == 2
        ylabel('Remote Insulin[1/min]')
    elseif i == 3
        ylabel('Subcutaneou Insulin(Isc1)[pmol/l]')
    elseif i == 4
        ylabel('Subcutaneou Insulin(Isc2)[pmol/l]')
    elseif i == 5
        ylabel('Plasma Insulin[pmol/l]')
    end
    % ax = gca;
    % saveFig(ax,'png')
end


close all
figure
a1 = bar(xlabel,obs,1);
legend({'4日間','4日目のみ'},'FontSize',15,'FontWeight','bold')
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;


figure
a2 = bar(xlabel,ekf,1);
legend({'4日間','4日目のみ'},'FontSize',15,'FontWeight','bold')
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
yticks(0:0.5:2)

figure
b1 = bar(xlabel,day,1);
legend({'Observer','EKF'},'FontSize',15,'FontWeight','bold')
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;

figure
b2 = bar(xlabel,day4,1);
legend({'Observer','EKF'},'FontSize',15,'FontWeight','bold')
h_axes = gca;
h_axes.XAxis.FontSize = 20;
h_axes.YAxis.FontSize = 20;
yticks(0:0.5:2)