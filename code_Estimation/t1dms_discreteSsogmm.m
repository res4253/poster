close all
clear
clc

% Path = ['./t1dms_data/:',...
%     './t1dms_data/subjects/:',...
%     './t1dms_data/simu_data/:'];

Path = [
    '.\t1dms_data\subjects\;',...
    '.\t1dms_data\simu_data\;'];
addpath(Path);


%load data file
load('temp.mat') %サンプリングデータのロード

load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])
p = p_esti;

%sampling time
N = size(data_ssogmm.ts,2);
N_step = 1:N;

%data
G = data.G.signals.values';

Xs = data_ssogmm.xs;
% Xs = x_esti;
%feedback
% Y = G(data_ssogmm.ts);
Y = data_ssogmm.Gs;

%%input
u_m = data_ssogmm.meal;
u_i = data_ssogmm.insulin;
modes = data_ssogmm.modes;

%Ra
R_a = get_Ra(Xs(4,:),p,modes,N);

U = [u_m;u_i;R_a;modes];


%%
Q = 1e-2; %システム雑音
R = 0.006; %観測雑音
% v = randn(1,N)*sqrtm(Q);
% w = randn(1,N)*sqrtm(R);


%%ssogmm
dic = dic_ssogmmCL(p,Ts,N,U);
x = dic.dicrete(Xs,"ZOH");


%kalman
x_ini = x_esti(:,500);
x_ini = x_ini([1 2 5 6 7]);


A = dic.A;
A_di = dic.A_di;
A_dx = dic.A_dx;
C = [1 0 0 0 0];

xhat = zeros(5,N);

a = [10,1e-6,1000,1000,500];
P = diag(a);
xhat(:,1) = x_ini;

B = [0.1, 1e-6, 100, 100, 25];
B_diag = diag(B);
for k = 2:N

     xhatm = dic.dynamics(xhat(:,k-1),k);
     
     A_tmp = [zeros(3,2) A_di];
     A_temp =[A(xhat(1,k-1),xhat(2,k-1)) ; A_tmp];
     % Pm = A(xhat(1,k-1),xhat(2,k-1))*P*A(xhat(1,k-1),xhat(2,k-1)) + b*Q*b';
     Pm = A_temp*P*A_temp' + B_diag*Q*B_diag';

     g = Pm*C'/(C*Pm*C' + R);

     xhat(:,k) = xhatm + g*(Y(:,k) - C*xhatm);
     P = (eye(5) - g*C)*Pm;
end


% f=figure;
% plot(Y,'LineWidth',2)
% hold on
% plot(x(1,:),'.')
% plot(xhat(1,:),'.')
% legend('uva/padova','x','xhat')


ssogmm = data_ssogmm.xs([1 2 5 6 7],:);
for i=1:5
    figure
    if i==1
    % plot(G)  
    hold on
    plot(data_ssogmm.ts,xhat(i,:),'.')
    plot(data_ssogmm.ts,ssogmm(i,:))
    legend('xhat','x')
    else
    hold on
    plot(data_ssogmm.ts,xhat(i,:),'.')
    plot(data_ssogmm.ts,ssogmm(i,:))
    legend('xhat','x')
    end
end


d = floor(N - N/4);
rmse(xhat(1,d:end)',Y(1,d:end)')

for i = 1:5
    figure(i)
    % hold on
    % plot(data_ssogmm.ts,ssogmm5(i,:))
    % legend
    grid on
    % ユーザーにファイル名を入力してもらう
    % filename = input('保存するファイル名を入力してください（拡張子は不要）: ', 's');
    if i == 1
        filename = 'g';
    elseif i == 2
        filename = 'x';
    elseif i == 3
        filename = 'i1';
    elseif i == 4
        filename = 'i2';
    elseif i == 5
        filename = 'ip';
    end

    % fig ファイルとして保存
    % savefig([filename, '.png']);
    ax = gca;
    saveFig(ax,'png','name',filename)
end