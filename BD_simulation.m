clear all; close all; clc;format compact;format shortG;
name = 'ryusiwon';
%name = 'hannu';
%cd(['C:\Users\',name,'\Dropbox\1_Study\2. TA\2019-1 Topics\TA session material\SImulation\BD']);

%A = importdata("c:\data_temp\M4_E_nonmissing.csv")
A = csvread('c:\data_temp\M4_E_nonmissing.csv',1)
v6 = A(:,6);
v3= A(:,3);
x1 = v3(v6==1);
x2 = v3(v6==0);
%%
% Set parameters
global R B ngrid alpha n1 n2
alpha   = 0.05;
R          = 1; % number of simulation
B          = 100; % number of bootstrap
ngrid   = 50;
n1      = 100;
n2      = 100;

% Set design [mean of x1, sd of x1, mean of x2, sd of x2]
% H0 : x1 weakly sd x2 <=> F1 <= F2
% BDstat = sqrt(-)*sup(F1-F2)
design1 = [0 1 0 1]; %FSD, SSD
design2 = [-0.1 1 0 1] ; %FSD
design3 = [0.5 1 0 1]; % SSD
design4 = [0 0.5 0 1]; %None
design = [design1;design2;design3;design4];  

% Graphing
% for m = 1 : 4
%     grid = linspace(-4,4,ngrid)';
%     
%     syms x_sym
%     cdfn = @(x,m,s) cdf('normal',x,m,s);
%     cdf1 = @(x) cdfn(x,design(m,1),design(m,2));
%     cdf2 = @(x) cdfn(x,design(m,3),design(m,4));
%     icdf1 = matlabFunction(int(cdf1(x_sym)));
%     icdf2 = matlabFunction(int(cdf2(x_sym)));
% 
%     Y1 = [cdf1(grid),cdf2(grid)];
%     Y2 = [icdf1(grid),icdf2(grid)];
%     
%     fig1 = figure(1)
%     subplot(2,2,m)
%     hold on
%     plot1 = plot(grid,Y1,'Color',[0 0 0]);
%             name1 = ['X1~N(',num2str(design(m,1)),',',num2str(design(m,2)),')'];
%             name2 = ['X2~N(',num2str(design(m,3)),',',num2str(design(m,4)),')'];
%             set(plot1(1),'DisplayName',name1,'LineStyle','--');
%             set(plot1(2),'DisplayName',name2);
%     legend1 = legend('show','Location','southeast');
%     hold off
%     set(fig1,'InnerPosition',[850,650,800,600])
%     
%     fig2 = figure(2)
%     subplot(2,2,m)
%     hold on
%     plot1 = plot(grid,Y2,'Color',[0 0 0]);
%             name1 = ['X1~N(',num2str(design(m,1)),',',num2str(design(m,2)),')'];
%             name2 = ['X2~N(',num2str(design(m,3)),',',num2str(design(m,4)),')'];
%             set(plot1(1),'DisplayName',name1,'LineStyle','--');
%             set(plot1(2),'DisplayName',name2);
%     legend1 = legend('show','Location','southeast');
%     hold off
%     set(fig2,'InnerPosition',[850,650,800,600])
% end


%Test
m = 1
sample1 = design(m,1)+design(m,2)*randn(n1,1,1,1);
sample2 = design(m,3)+design(m,4)*randn(n2,1,1,1);
grid = linspace(min(min([sample1;sample2])), max(max([sample1;sample2])),ngrid);
SDorder = 1;
[bd,pval]  = bdtest(sample1,sample2,grid,SDorder)


%% Simulation
tic;
R = 100; % Number of simulation
clear bd_R  pval_R
for SDorder = [1]
    for m = [1]
        sample1 = design(m,1)+design(m,2)*randn(n1,1,1,R);
        sample2 = design(m,3)+design(m,4)*randn(n2,1,1,R);
        grid = linspace(min(min([sample1;sample2])), max(max([sample1;sample2])),ngrid);
        [bd_R(:,m,SDorder),pval_R(:,m,SDorder,:)]  = bdtest(sample1,sample2,grid,SDorder);
    end
    disp(['Design ',num2str(m),' of SD order ',num2str(SDorder), ' Done'])
end
if R > 1
    rej_prob = mean(pval_R<alpha,4);
    disp(['Rejection probabilty = '])
    disp(rej_prob)
else
    disp('bd = ') 
    disp(bd_R)
    disp('pvalues = ')
    disp(pval_R)
end
save('result_BD_new')
toc;


function [bd,pval] = bdtest(sample1,sample2,grid,SDorder)
global R B
n1 = size(sample1,1); n2 = size(sample2,1);

operator = @(X,z) (X<=z).*(z-X).^(SDorder-1)/factorial(SDorder-1); %integral operator
ecdf    = @(X,z) mean(bsxfun(operator,X,z));
stat    = @(z) sqrt(n1*n2/(n1+n2))*max(ecdf(sample1,z)-ecdf(sample2,z),[],2);
bd      = stat(grid);

% Multiplier method
% multiplier1 : BD(2003), (6)
% multiplier2 : BD(2003), (5)
lambda = n2/(n1+n2);

op_J    = @(X,z) sqrt(size(X,1))*mean([operator(X,z) - ecdf(X,z)].*randn(size(X,1),1,B,R),1);
multiplier1 = @(z) max(sqrt(lambda)*op_J(sample1,z) - sqrt(1-lambda)*op_J(sample2,z),[],2);
multiplier2 = @(z) max(op_J(sample2,z),[],2);

% 2. Bootstrap
% 1) Separating
b1sample1 = sample1(randi([1,n1],n1,1,B,R));
b1sample2 = sample2(randi([1,n2],n2,1,B,R));

% 2) Pooling
sample_p    = [sample1; sample2];
b2sample1        = sample_p(randi([1,n1+n2],n1,1,B,R)) ;
b2sample2        = sample_p(randi([1,n1+n2],n2,1,B,R));

% 3) Statistics
% boot1 : using sample2(BD(2003), (9))
% boot2 : pooled bootstrap(BD(2003), (10))
% boot3 : recentered bootstrap(BD(2003), (11))
boot1 = @(z) sqrt(n2)*max(ecdf(b1sample2,z) - ecdf(sample2,z),[],2);
boot2 = @(z) sqrt(n1*n2/(n1+n2))*max(ecdf(b2sample1,z) - ecdf(b2sample2,z),[],2);
boot3 = @(z) sqrt(n1*n2/(n1+n2))*max((ecdf(b1sample1,z)-ecdf(sample1,z)) -(ecdf(b1sample2,z)-ecdf(sample2,z)),[],2);

% 2. Results
pval = [
    mean(multiplier1(grid) > bd,3);
    mean(multiplier2(grid) > bd,3);
    mean(boot1(grid) > bd,3);
    mean(boot2(grid) > bd,3);
    mean(boot3(grid) > bd,3)
    ];
end




