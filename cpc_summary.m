function [] = cpc_summary(ID)
global dataPath Msamples burn Nchains datachoice filename multipleconst
global simu MinNtimesw  inc_int rep model estAlpha inc_ints

id = str2double(num2str(ID)); % 1-90
rng('default'); rng(id*241);
rep = ceil(id/3);
model = id - (rep-1)*3;
datachoice = 1;
simu = 1;
inc_int = 0; inc_ints = 1;
if simu == 1
    Msamples = 1e3; burn = 0; Nchains = 3;
else
    Msamples = 1e3; burn = 0; Nchains = 4;
end
estAlpha = 1;
multipleconst = 1;
dataPath = './';
if datachoice == 1
    MinNtimesw = 7;
elseif datachoice == 3
    MinNtimesw = 3;
end
filename = strcat(dataPath,'figuremodel',num2str(model),'_',num2str(rep),'.ps');
ConvergenceDiagnostics(model, rep)
checkmodelpro(model, rep)
end

function [] = ConvergenceDiagnostics(model, rep)
%function to check the convergence of MC samples
global dataPath Msamples burn Nchains filename datachoice inc_int simu timesize
if simu == 1
    load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat'))
    timesize = size(data,1)-1;
end
if datachoice == 1
    %set monitoring sites
    monitorsite=[1 12 18 27 37 46];
    %set monitoring times
    monitortime=[1 14 24 35];
elseif datachoice == 2
    monitorsite=[57 12 87];
    monitortime=[8 20 35];
elseif datachoice == 3
    monitorsite=[12 16 24 52 61 67];
    monitortime=[1 4 8 12 17];
end

nsite = numel(monitorsite); ntime = numel(monitortime);

%D is number of parameters
%M is number of MCMC samples
%N is number of chains
%NxDxM matrix to put in psrf.m function

stlen = numel(monitorsite)*numel(monitortime);
Dparameters = 2+nsite*2+stlen;
if  inc_int == 1
    Dparameters = Dparameters + numel(monitorsite)*numel(monitortime);
end

MCSamples = zeros(Nchains,Dparameters+2,Msamples);

tmp = cell(Msamples,Nchains);
ds = zeros(Msamples+burn,Nchains);
n2loglik = zeros(Msamples+burn,Nchains);

for inits = 1:Nchains
    load(strcat(dataPath,'MCMCoutputCh',num2str(model),'_ini',num2str(inits),'_',num2str(rep),'.mat'))
    for iter = 1:length(xall) %(length(xall) - Msamples + 1)
        if iter > burn
            tmp{iter-burn, inits} = xall{iter};
        end
        Nstar = length(xall{iter});
        ds(iter, inits) = Nstar;
        loglik = 0;
        for r = 1:Nstar
            nr = length(xall{iter}{r}.Siteincluster);
            x = xall{iter}{r};
            obj = computefromInt(x.NpointinTimeIn);
            for s = 1:nr
                obsY = data(2:end, x.Siteincluster(s));
                tmpDelta = obsY - obj.X1*x.parameter' - obj.X0*x.beta0{s}';
                loglik = loglik - 0.5*length(obsY)*log(2*pi*x.sigma2(s)) - 0.5/x.sigma2(s)*(sum(tmpDelta.^2));
            end
        end
        n2loglik(iter, inits) = -2*loglik;
    end
end
xall = tmp; clear tmp;

plot(n2loglik)
title('-2loglikelihood statistic for monitoring convergence');
legend('chain 1','chain 2','chain 3','Location','Best')
orient landscape
print('-painters', '-dpsc2', '-r600', filename)

for iter = 1:Msamples
    for nch = 1:Nchains
        MCSamples(nch,1:2,iter) = [xall{iter,nch}{1,1}.tau2, xall{iter,nch}{1,1}.gamma];
        MCSamples(nch,Dparameters+(1:2),iter) = [xall{iter,nch}{1,1}.lambda, xall{iter,nch}{1,1}.alpha];
        Nstar = size(xall{iter,nch},1);% define N*: the number of cluster
        for st = 1:nsite
            onesite = monitorsite(st);
            for ncl = 1:Nstar
                marks = abs(xall{iter,nch}{ncl,1}.Siteincluster-onesite);
                if ~min(marks)
                    sind = find(~marks);
                    MCSamples(nch,2+st,iter) = xall{iter,nch}{ncl,1}.sigma2(sind);  %sigma2
                    MCSamples(nch,2+st+nsite,iter) = xall{iter,nch}{ncl,1}.u(sind); %u
                    for tm = 1:ntime
                        indx = tm+(st-1)*ntime;
                        onetime = monitortime(tm);
                        TmpTimeInt = [0 xall{iter,nch}{ncl,1}.NpointinTimeIn];
                        SizeInt = numel(TmpTimeInt);
                        for titer=2:SizeInt
                            if and(sum(TmpTimeInt(1:(titer-1)))<onetime,sum(TmpTimeInt(1:titer))>=onetime)
                                timeindx=titer-1;
                            end
                        end
                        MCSamples(nch,2+nsite*2+indx,iter) = xall{iter,nch}{ncl,1}.parameter(timeindx);
                        if inc_int == 1
                            MCSamples(nch,2+nsite*2+stlen+indx,iter) = xall{iter,nch}{ncl,1}.beta0(timeindx);
                        end
                    end
                end
            end
        end
    end
end

R = psrf(MCSamples);
display(R)

plot(reshape(MCSamples(:,1,:),[Nchains,Msamples])')
series = reshape(MCSamples(:,1,:),[Nchains*Msamples,1]);
[l,u] = FindHPDset(series,0.95,[]);
title({strcat('\tau^2 with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
legend('chain 1','chain 2','chain 3','chain 4','Location','Best')
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)
hist(series); title('histogram of \tau^2 for pooled chains');
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

plot(reshape(MCSamples(:,2,:),[Nchains,Msamples])')
series = reshape(MCSamples(:,2,:),[Nchains*Msamples 1]);
[l,u] = FindHPDset(series,0.95,[]);
title({strcat('\gamma with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)
hist(series); title('histogram of \gamma for pooled chains');
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

plot(reshape(MCSamples(:,Dparameters+1,:),[Nchains,Msamples])')
series = reshape(MCSamples(:,Dparameters+1,:),[Nchains*Msamples,1]);
[l,u] = FindHPDset(series,0.95,[]);
title({strcat('\lambda with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)
hist(series); title('histogram of \lambda for pooled chains');
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

plot(reshape(MCSamples(:,Dparameters+2,:),[Nchains,Msamples])')
series = reshape(MCSamples(:,Dparameters+2,:),[Nchains*Msamples,1]);
[l,u] = FindHPDset(series,0.95,[]);
title({strcat('\alpha_0 with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)
hist(series); title('histogram of \alpha_0 for pooled chains');
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

for j = 1:nsite
    plot(reshape(MCSamples(:,2+j,:),  [Nchains Msamples])')
    series = reshape(MCSamples(:,2+j,:),[Nchains*Msamples 1]);
    [l,u] = FindHPDset(series,0.95,[]);
    title({strcat('\sigma^2 for site ',num2str(monitorsite(j)),' with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
    orient landscape
    print('-painters', '-dpsc2', '-r600', '-append', filename)
end

for j = 1:nsite
    plot(reshape(MCSamples(:,2+j+nsite,:),  [Nchains Msamples])')
    series = reshape(MCSamples(:,2+j+nsite,:),[Nchains*Msamples 1]);
    [l,u] = FindHPDset(series,0.95,[]);
    title({strcat('u for site ',num2str(monitorsite(j)),' with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
    orient landscape
    print('-painters', '-dpsc2', '-r600', '-append', filename)
end

for j = 1:nsite
    for k = 1:ntime
        indx = k+(j-1)*ntime;
        plot(reshape(MCSamples(:,2+nsite*2+indx,:),  [Nchains Msamples])')
        series = reshape(MCSamples(:,2+nsite*2+indx,:),[Nchains*Msamples 1]);
        [l,u] = FindHPDset(series,0.95,[]);
        title({strcat('\beta for site ',num2str(monitorsite(j)),' time ',num2str(monitortime(k)),'with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
        orient landscape
        print('-painters', '-dpsc2', '-r600', '-append', filename)
    end
end

if inc_int
    for j = 1:nsite
        for k = 1:ntime
            indx = k+(j-1)*ntime;
            plot(reshape(MCSamples(:,2+nsite*2+stlen+indx,:),  [Nchains Msamples])')
            series = reshape(MCSamples(:,2+nsite*2+stlen+indx,:),[Nchains*Msamples 1]);
            [l,u] = FindHPDset(series,0.95,[]);
            title({strcat('intercept for site ',num2str(monitorsite(j)),' time ',num2str(monitortime(k)),'with mean, median, 95% HPD set:  ',num2str([mean(series), median(series)]));strcat('[',num2str([l u]),']')});
            orient landscape
            print('-painters', '-dpsc2', '-r600', '-append', filename)
        end
    end
end

save(strcat(dataPath,'MCMCoutputModel',num2str(model),'_',num2str(rep),'.mat'),'xall', 'R')
display('done')
end

function [] = checkmodelpro(model, rep)
global t timesize simu inc_int MinNtimesw inc_ints
global sitesize dataPath Msamples Nchains datachoice filename true_bvec Clustergroup0 true_cpt

if datachoice == 1
    load('neighbor.txt')
else
    load('iowa_neighbor.mat'); neighbor = NB; %larger
end

if simu == 1
    load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat'))   % get true_bvec
    data = data(2:end,:); %#ok<NODEF>
else
    if datachoice == 1
        load('data.txt') % for real data
        data = log(data(2:end,:)); %#ok<NODEF>
        u0 = zeros(1,size(data,2));
        if inc_ints ~= 1
            for j = 1:size(data,2)
                ay = data(1:MinNtimesw,j); ax = 1:MinNtimesw;
                betas =regress(ay, [ones(length(ax),1), ax']);
                u0(j) = betas(1);   % mean vector of random effects
            end
            data = data - repmat(u0, [size(data,1) 1]);
        end
        
    elseif  datachoice == 2
        data = load('iowa09-12.txt'); data = data'; data = data(5:35,:);
    elseif  datachoice == 3
        data = load('yeast cell data.txt');
    end
end

t = cat(1,1:size(data,2),data);
if inc_int ~= 1 && inc_ints ~= 1
    GenerateRSSandbetamodel1(t, model, rep);
end

data = t(2:end,:);
site = t(1,:);
sitesize = size(site,2);
% timesize = size(data,1);

HPDp = 0.95;
pp = 0.05;
load(strcat(dataPath,'MCMCoutputModel',num2str(model),'_',num2str(rep),'.mat'))
for j = 1:Msamples
    for chs = 1:Nchains
        xallIf{j+Msamples*(chs-1),1} = xall{j,chs};
    end
end
clear xall

ClusteIndex = zeros(sitesize,1);
TClustersize = size(xallIf,1);
for i = 1:TClustersize,
    ClusteIndex(size(xallIf{i,1},1)) = ClusteIndex(size(xallIf{i,1},1))+1;
end
[val,c] = max(ClusteIndex);
afteriteration=1;
totIter = size(xallIf,1);

tau2 = zeros(1, totIter); gamma = zeros(1,totIter);
lambda = zeros(1, totIter); alpha = zeros(1,totIter);
for loop = 1:totIter
    tau2(loop) = xallIf{loop,1}{1}.tau2;
    gamma(loop) = xallIf{loop,1}{1}.gamma;
    lambda(loop) = xallIf{loop,1}{1}.lambda;
    alpha(loop) = xallIf{loop,1}{1}.alpha;
end
[l,u] = FindHPDset(gamma',HPDp,[]); HPD.gamma = [l mean(gamma) median(gamma) u];
[l,u] = FindHPDset(tau2',HPDp,[]); HPD.tau2 = [l mean(tau2) median(tau2) u];
[l,u] = FindHPDset(lambda',HPDp,[]); HPD.lambda = [l mean(lambda) median(lambda) u];
[l,u] = FindHPDset(alpha',HPDp,[]); HPD.alpha = [l mean(alpha) median(alpha) u];

l = quantile(gamma,(1-HPDp)/2); u = quantile(gamma,1-(1-HPDp)/2);
HPD.gammaCI = [l mean(gamma) u];
l = quantile(tau2,(1-HPDp)/2); u = quantile(tau2,1-(1-HPDp)/2);
HPD.tau2CI = [l mean(tau2) u];
l = quantile(lambda,(1-HPDp)/2); u = quantile(lambda,1-(1-HPDp)/2);
HPD.lambdaCI = [l mean(lambda) u];
l = quantile(alpha,(1-HPDp)/2); u = quantile(alpha,1-(1-HPDp)/2);
HPD.alphaCI = [l mean(alpha) u];

% cluster information
[Clustergroup,SS,SS2] = makelinkage(xallIf,c);
if simu == 1 % create the 2 by 2 cross-classification table and measures
    clusterMeasure = computefromTable(Clustergroup, Clustergroup0);
    display('cross-classification results:')
    display(clusterMeasure.tab)
    display(clusterMeasure)
    outH = getDistance2(data, Clustergroup0, Clustergroup);
    scatter(clusterMeasure.sensitivity, clusterMeasure.specificity, 60,'r','filled');
    hold off;
else
    hold on; outH = getDistance2(data, [], Clustergroup); hold off;
end
if simu == 0
    save('outH.mat','outH')
else
    save(strcat('outH_',num2str(rep),'.mat'),'outH')
end
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

PsOutput = cell(c,1);
fitMeasure.Pred = zeros(timesize, sitesize);
fitMeasure.Drep = zeros(totIter-afteriteration+1, sitesize);
fitMeasure.Dhat = zeros(1, sitesize);
fitMeasure.SSE = zeros(1, sitesize);
fitMeasure.MSE = zeros(1, sitesize);
fitMeasure.BIC = zeros(1, sitesize);
fitMeasure.DIC = zeros(1, sitesize);

cptMeasure.distcpt = zeros(timesize,sitesize);  % mean posterior distribution of changepoint at each time, each site
cptMeasure.HPD = cell(1,sitesize);
cptMeasure.HPDp = cell(1,sitesize);
cptMeasure.KL = zeros(1,sitesize); % KL divergence between posterior distribution of changepoint and the true mass
cptMeasure.entropy = zeros(1,sitesize); % entropy of posterior
cptMeasure.relentropy = zeros(1,sitesize); % relative entropy of posterior
cptMeasure.accuracy = zeros(1,sitesize); % mean accuracy rate for detection of true changepoints
cptMeasure.betadiff = zeros(1,sitesize); % difference of beta
cptMeasure.betabwHPD = zeros(1,sitesize); % HPD band width of difference
cptMeasure.betabwCI = zeros(1,sitesize); % CI band width of difference

% randn('state',267);
Rdata = cell(1,c);
for i = 1:c,
    monitorsite = Clustergroup{i,1};
    [PsOutput{i,1},fitMeasure,cptMeasure,Rdata{i}] = getPSOutput(model, xallIf,afteriteration,t,monitorsite,totIter,HPDp,pp,...
        fitMeasure, cptMeasure);
    orient landscape
    %     print('-painters', '-depsc2', '-r600', strcat(num2str(i),filename))
    print('-painters', '-dpsc2', '-r600', '-append', filename)
end
fitMeasure.Dave = mean(-2*fitMeasure.Drep, 1); fitMeasure.Dhat = -2*fitMeasure.Dhat;
fitMeasure.DIC = 2*fitMeasure.Dave - fitMeasure.Dhat;
display('DIC:')
display(fitMeasure.DIC)
fitMeasure.Drep = []; % no need to store

% Name = char('Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Washington DC','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky',...
% 'Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire	',...
% 'New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina',...
% 'South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming');

%%% probability of number of clusters
postCluster = [(1:11)' ClusteIndex(1:11)/(Msamples * Nchains)];
display('posterior distribution of clusters:')
display(postCluster)
tmpvec = postCluster(postCluster(:,2) ~=0, 2);
clusterMeasure.SS = SS;
clusterMeasure.SS2 = SS2;
clusterMeasure.postCluster = postCluster;
clusterMeasure.entropy = - sum(tmpvec.*log(tmpvec));
if length(tmpvec) == 1
    clusterMeasure.relentropy = 1;
else
    clusterMeasure.relentropy = 1 - clusterMeasure.entropy/log(length(tmpvec));
end
display('entropy:')
display(clusterMeasure.entropy)
save(strcat(dataPath,'Rdata',num2str(model),'_',num2str(rep),'.mat'),'Rdata')
save(strcat(dataPath,'PsOutputModel',num2str(model),'_',num2str(rep),'.mat'),...
    'PsOutput','fitMeasure','HPD','R','clusterMeasure','cptMeasure','Clustergroup','outH')
close all
end

function [psoutput,fitMeasure,cptMeasure,Rdata] = ...
    getPSOutput(model, xalliter,afteriter,alldata,allmonitorsite,totIter,HPDp,pp,fitMeasure,cptMeasure)
global timesize dataPath  datachoice rep sitesize  MinNtimesw simu
global multipleconst inc_int inc_ints true_bvec  Clustergroup0 estAlpha true_cpt
if inc_int ~= 1 && inc_ints ~= 1
    load(strcat(dataPath,'Ressandbeta',num2str(model),'_',num2str(rep),'.mat'))  % beta0s
end
monitorparameterindex=[1]; thetatemp=[];
beta1temp=cell(size(allmonitorsite,2),1); beta0temp=cell(size(allmonitorsite,2),1);
sigmavaltemp=cell(size(allmonitorsite,2),1); uvaltemp=cell(size(allmonitorsite,2),1);  c=[];
Name = cell(1, sitesize);
for inds = 1:sitesize
    Name{inds} = num2str(inds);
end
if datachoice == 1
    Name=char('Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Washington DC','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky',...
        'Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire	',...
        'New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina',...
        'South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming');
end
figure; hold on;
miny=((1./multipleconst).*(min(min(alldata(2:end,allmonitorsite)))-5));
maxy=((1./multipleconst).*(max(max(alldata(2:end,allmonitorsite)))+5));
psoutput.sigma=zeros(size(allmonitorsite,2),1);
psoutput.theta=cell(size(allmonitorsite,2),1);
psoutput.changetime=cell(size(allmonitorsite,2),1);
psoutput.NTimeInterals=zeros(size(allmonitorsite,2),1);
psoutput.NpointinTimeIn=cell(size(allmonitorsite,2),1);
psoutput.Siteincluster=allmonitorsite;
psoutput.HPDsetdata=cell(size(allmonitorsite,2),1);
psoutput.yvalue=cell(size(allmonitorsite,2),1);
psoutput.sepredT=cell(size(allmonitorsite,2),1);
psoutput.betadifference=cell(size(allmonitorsite,2),1);
psoutput.betadifferenceHPD=cell(size(allmonitorsite,2),2);
psoutput.betadifferenceCI=cell(size(allmonitorsite,2),2);
psoutput.numTCounts=cell(size(allmonitorsite,2),1);
psoutput.numTintervalCounts=cell(size(allmonitorsite,2),1);
psoutput.thetaCI = cell(size(allmonitorsite,2),2);

yy = zeros(totIter-afteriter+1, size(allmonitorsite,2));
nr = size(allmonitorsite,2);
Rdata.sites = allmonitorsite;
Rdata.ObsY = zeros(timesize,nr);
Rdata.miny = zeros(1,nr);
Rdata.maxy = zeros(1,nr);
% Rdata.Title = zeros(1,nr);
Rdata.sigma2 = zeros(1,nr);
Rdata.yvalueupperCI = zeros(timesize,nr);
Rdata.yvaluelowerCI = zeros(timesize,nr);
Rdata.cpt = cell(1,nr);

for monidex = 1:size(allmonitorsite,2)
    if length(allmonitorsite) == 8
        if monidex ~= 8
            subplot(3,4,monidex);
        else
            subplot(3,4,monidex+1);
        end
    elseif length(allmonitorsite) == 16
        subplot(4,4,monidex);
    else
        subplot(floor(sqrt(size(allmonitorsite,2)))+1,floor(sqrt(size(allmonitorsite,2)))+1,monidex);
    end
    obsY = ((1./multipleconst).*alldata(2:end,allmonitorsite(monidex)));
    plot(1:timesize,obsY,...
        'o','MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
    hold on;
    axis([0 timesize+1 miny maxy]);
    Rdata.ObsY(:,monidex) = obsY;
    if datachoice ~= 1
        title(Name{allmonitorsite(monidex)}, 'Fontsize', 16)
    else
        title(Name(allmonitorsite(monidex),:), 'Fontsize', 16)
    end
    % Rdata.Title(monidex) = Name(allmonitorsite(monidex),:);
    AllHPDsetdata = zeros(timesize,(totIter-afteriter+1));
    totalTinterval = [(1:1:timesize)',zeros(timesize,1)];
    carddiff = 0;
    for ss = afteriter:totIter,
        Nstar = size(xalliter{ss,:},1);% define N*: the number of cluster
        HPDsetdata=[];
        for j = 1:Nstar,
            if (min(abs(xalliter{ss,1}{j,1}.Siteincluster-allmonitorsite(monidex)))==0)
                % allmonitorsite(monidex) lies in jth cluster at ssth iteration
                c(ss,monidex) = j;
                tmp0 = cumsum([1 xalliter{ss,1}{j,1}.NpointinTimeIn]);
                obs_cpt = tmp0(2:(length(tmp0)-1));  % observed changepoint
                if simu == 1
                    for i0 = 1:size(Clustergroup0)
                        if any(Clustergroup0{i0} == allmonitorsite(monidex))
                            k0 = i0;
                        end
                    end
                    if length(obs_cpt)>length(true_cpt{k0})
                        carddiff = carddiff + length(obs_cpt)-length(true_cpt{k0});
                    end
                end
                cptMeasure.distcpt(obs_cpt, allmonitorsite(monidex)) = cptMeasure.distcpt(obs_cpt, allmonitorsite(monidex)) + 1;
                [xx,yy(ss,monidex)] = min(abs(xalliter{ss,1}{j,1}.Siteincluster-allmonitorsite(monidex)));
                % yy: the position  of the site allmonitorsite(monidex) in the cluster j
                TimeIntervalHPD = [0 cumsum(xalliter{ss,1}{j,1}.NpointinTimeIn,2)];
                BHPD = Combination(xalliter{ss,1}{j,1}.NTimeIntervals,timesize-xalliter{ss,1}{j,1}.NTimeIntervals*MinNtimesw);
                [xxxtmep,yyytemp] = min(sum((BHPD-ones(size(BHPD,1),1)*xalliter{ss,1}{j,1}.NpointinTimeIn).^2,2));
                obj = computefromInt(xalliter{ss,1}{j,1}.NpointinTimeIn);
                tmpSigma2 = (1./multipleconst.^2).*xalliter{ss,1}{j,1}.sigma2(yy(ss,monidex));
                tmpU = (1./multipleconst).*xalliter{ss,1}{j,1}.u(yy(ss,monidex));
                tmpDelta = obsY - tmpU - ...
                    (1./multipleconst).*obj.X1*xalliter{ss,1}{j,1}.parameter';
                %                 Drep(ss) = Drep(ss) - 0.5*(timesize - length(TimeIntervalHPD) + 1)*log(2*pi*tmpSigma2)...
                %                     - 0.5*log(det(obj.X0' * obj.X0)) - 0.5/tmpSigma2*tmpDelta'*obj.Sigma*tmpDelta;
                if inc_int ~= 1 && inc_ints ~= 1
                    Betazeros = beta0s{xalliter{ss,1}{j,1}.NTimeIntervals,allmonitorsite(monidex)}(yyytemp,:);
                    if estAlpha == 1 && obj.n1 < timesize
                        Betazeros(2:end) = mvnrnd(inv(obj.X0' * obj.X0)*obj.X0'*tmpDelta((obj.n1+1):end), tmpSigma2*inv(obj.X0' * obj.X0));
                    end
                elseif inc_int == 1
                    Betazeros = xalliter{ss,1}{j,1}.beta0;
                elseif inc_ints == 1
                    Betazeros = xalliter{ss,1}{j,1}.beta0{yy(ss,monidex)};
                end
                if inc_ints ~= 1
                    fitMeasure.Drep(ss,allmonitorsite(monidex)) = fitMeasure.Drep(ss,allmonitorsite(monidex)) - 0.5*timesize*log(2*pi*tmpSigma2) - 0.5/tmpSigma2*(sum(tmpDelta(1:obj.n1).^2) + ...
                        sum((tmpDelta((obj.n1+1):end) - obj.X0*Betazeros(2:end)').^2));
                else
                    fitMeasure.Drep(ss,allmonitorsite(monidex)) = fitMeasure.Drep(ss,allmonitorsite(monidex)) - 0.5*timesize*log(2*pi*tmpSigma2) - 0.5/tmpSigma2*(sum(tmpDelta(1:obj.n1).^2) + ...
                        sum((tmpDelta - obj.X0*Betazeros').^2));
                end
                
                for i = 1:size(TimeIntervalHPD,2)-1,
                    HPDsettemp = [];
                    timeaxisHPD = [TimeIntervalHPD(i)+1:1:TimeIntervalHPD(i+1)];
                    tmpmean = timeaxisHPD.*xalliter{ss,1}{j,1}.parameter(i) + xalliter{ss,1}{j,1}.u(yy(ss,monidex)) + Betazeros(i);
                    
                    for tm = 1:numel(tmpmean)
                        HPDsettemp(tm)=(1./multipleconst).*(normrnd(tmpmean(tm),sqrt(xalliter{ss,1}{j,1}.sigma2(yy(ss,monidex)))));
                    end
                    HPDsetdata = [HPDsetdata HPDsettemp];
                end
                AllHPDsetdata(:,ss)= HPDsetdata';
                
                utemp(ss,monidex) = tmpU;
                sigmatemp(ss,monidex) = tmpSigma2;
                thetatemp(ss,monidex) = (1./multipleconst).*(xalliter{ss,1}{j,1}.parameter(monitorparameterindex));
                timeintervaltemp = size(xalliter{ss,1}{j,1}.NpointinTimeIn, 2); % # of time intervals
                totalTinterval(timeintervaltemp,2) = totalTinterval(timeintervaltemp,2) + 1;
            else
            end
        end
    end
    
    carddiff = carddiff/(totIter-afteriter+1);
    tmpdistcpt = cptMeasure.distcpt(:, allmonitorsite(monidex))/(totIter-afteriter+1);
    cptMeasure.distcpt(:, allmonitorsite(monidex)) = tmpdistcpt/sum(tmpdistcpt);  % make it a probability
    [tmpdistcptS, I] = sort(cptMeasure.distcpt(:, allmonitorsite(monidex)), 'descend');
    ind = find(cumsum(tmpdistcptS) >= HPDp); ind = ind(1);
    cptMeasure.HPD{allmonitorsite(monidex)} = I(1:ind);
    cptMeasure.HPDp{allmonitorsite(monidex)} = cumsum(tmpdistcptS(1:ind));
    distQ = cptMeasure.distcpt(:, allmonitorsite(monidex)); distQ = distQ(distQ>0);
    cptMeasure.entropy(allmonitorsite(monidex)) = - sum(distQ.*log(distQ));
    if length(distQ) == 1
        cptMeasure.relentropy(allmonitorsite(monidex)) = 1;
    else
        cptMeasure.relentropy(allmonitorsite(monidex)) = 1 + sum(distQ.*log(distQ))/log(length(distQ));
    end
    if simu == 1
        distP0 = ones(1,length(true_cpt{k0}))/length(true_cpt{k0}); distQ0 = cptMeasure.distcpt(true_cpt{k0}, allmonitorsite(monidex)) + 1e-50;
        cptMeasure.KL(allmonitorsite(monidex)) = sum(distP0.*log(distP0./distQ0'));
        cptMeasure.accuracy(allmonitorsite(monidex)) = mean(tmpdistcpt(true_cpt{k0})) - carddiff/length(true_cpt{k0}); %+tmpdistcpt(true_cpt{k0}+1));
    end
    
    yvalueupperHPD = zeros(timesize,1);
    yvaluelowerHPD = zeros(timesize,1);
    yvaluelowerCI = zeros(timesize,1);
    yvalueupperCI = zeros(timesize,1);
    
    tmpYs2 = [];
    for sss = 1:timesize,
        [a,b] = FindHPDset(AllHPDsetdata(sss,:)',HPDp,[]);
        yvaluelowerCI(sss) = quantile(AllHPDsetdata(sss,:)',(1-HPDp)/2);
        yvalueupperCI(sss) = quantile(AllHPDsetdata(sss,:)',1-(1-HPDp)/2);
        if (size(a,2)==1)
            yvaluelowerHPD(sss) = a(1);
            yvalueupperHPD(sss) = b(1);
        else
            tempLength = zeros(size(a,2),1);
            for hindx = 1:size(a,2)
                tempLength(hindx) = b(hindx)-a(hindx);
            end
            [maxv,maxindx] = max(tempLength);
            yvaluelowerHPD(sss) = a(maxindx);
            yvalueupperHPD(sss) = b(maxindx);
        end
        %         for hpdi=1:size(a,2)
        %             plot(sss,(a(hpdi)),'^k'); hold on;
        %             plot(sss,(b(hpdi)),'vk'); hold on;
        %             tmpYs2 = [tmpYs2, a(hpdi), b(hpdi)];
        %         end
        tmpYs2 = [tmpYs2, yvaluelowerCI(sss), yvalueupperCI(sss)];
    end
    
    [test,maxtimeIn] = max(totalTinterval(:,2));
    psoutput.numTCounts{monidex} = totalTinterval; %total counts of MCsamples for each number of intervals
    psoutput.NTimeInterals(monidex) = maxtimeIn;
    B = Combination(maxtimeIn,timesize-maxtimeIn*MinNtimesw);
    sigmatempvalue = cell(size(B,1),1);  utempvalue=cell(size(B,1),1);
    testIndex = zeros(size(B,1),1);
    Beta1parameter = cell(size(B,1),1);
    Beta0parameter = cell(size(B,1),1);
    
    for ss = afteriter:totIter,
        if (size(xalliter{ss,1}{c(ss,monidex),1}.NpointinTimeIn,2) == maxtimeIn),
            % for those iterations with posterior mode of # of changepoints only
            [xxxtmep,yyytemp] = min(sum((B-ones(size(B,1),1)*xalliter{ss,1}{c(ss,monidex),1}.NpointinTimeIn).^2,2));
            sigmatempvalue{yyytemp,1} = [sigmatempvalue{yyytemp,1}; sigmatemp(ss,monidex)];
            utempvalue{yyytemp,1} = [utempvalue{yyytemp,1}; (utemp(ss,monidex))];
            testIndex(yyytemp) = testIndex(yyytemp)+1;
            Beta1parameter{yyytemp,1} = [Beta1parameter{yyytemp,1}; (1./multipleconst).*(xalliter{ss,1}{c(ss,monidex),1}.parameter)];
            if inc_int ~= 1 && inc_ints ~= 1
                tmpBeta0parameter = (1./multipleconst).*(beta0s{maxtimeIn,allmonitorsite(monidex)}(yyytemp,:));
                if estAlpha == 1
                    obj = computefromInt(xalliter{ss,1}{c(ss,monidex),1}.NpointinTimeIn);
                    if obj.n1 < timesize
                        tmpDelta = obsY - utemp(ss,monidex) - ...
                            (1./multipleconst).*obj.X1*xalliter{ss,1}{c(ss,monidex),1}.parameter';
                        tmpBeta0parameter(2:end) = mvnrnd(inv(obj.X0' * obj.X0)*obj.X0'*tmpDelta((obj.n1+1):end), sigmatemp(ss,monidex)*inv(obj.X0' * obj.X0));
                    end
                end
                Beta0parameter{yyytemp,1}=[Beta0parameter{yyytemp,1}; tmpBeta0parameter];
            elseif inc_int == 1
                Beta0parameter{yyytemp,1}=[Beta0parameter{yyytemp,1}; (1./multipleconst).*(xalliter{ss,1}{c(ss,monidex),1}.beta0)];
            elseif inc_ints == 1
                Beta0parameter{yyytemp,1}=[Beta0parameter{yyytemp,1}; (1./multipleconst).*(xalliter{ss,1}{c(ss,monidex),1}.beta0{yy(ss,monidex)})];
            end
        else
        end
    end
    
    [val,ComIndex] = max(testIndex);
    psoutput.numTintervalCounts{monidex} = [B testIndex];
    psoutput.NpointinTimeIn{monidex,1} = B(ComIndex,:);
    TimeInterval = [0 cumsum(B(ComIndex,:))];
    %v1temp{monidex,1}=vls{maxtimeIn,allmonitorsite(monidex)}(ComIndex,:);
    sigmavaltemp{monidex,1} = mean(sigmatempvalue{ComIndex,1});
    uvaltemp{monidex,1} = mean(utempvalue{ComIndex,1});
    beta1temp{monidex,1} = mean(Beta1parameter{ComIndex,1},1);
    
    beta0temp{monidex,1} = mean(Beta0parameter{ComIndex,1},1);
    psoutput.sigma2(monidex) = sigmavaltemp{monidex,1};
    psoutput.u(monidex) = uvaltemp{monidex,1};
    psoutput.theta{monidex,1} = beta1temp{monidex,1};
    
    obj = computefromInt(B(ComIndex,:));
    tmpSigma2 = mean(sigmatempvalue{ComIndex,1});
    tmpDelta = obsY - mean(utempvalue{ComIndex,1}) - obj.X1*mean(Beta1parameter{ComIndex,1},1)';
    %     Dhat = Dhat - 0.5*(timesize - length(B(ComIndex,:)) + 1)*log(2*pi*tmpSigma2)...
    %         - 0.5*log(det(obj.X0' * obj.X0)) - 0.5/tmpSigma2*tmpDelta'*obj.Sigma*tmpDelta;
    
    if inc_ints ~= 1
        fitMeasure.Dhat(allmonitorsite(monidex)) = fitMeasure.Dhat(allmonitorsite(monidex)) - 0.5*timesize*log(2*pi*tmpSigma2) - 0.5/tmpSigma2*(sum(tmpDelta(1:obj.n1).^2) + ...
            sum((tmpDelta((obj.n1+1):end) - obj.X0*beta0temp{monidex,1}(2:end)').^2));
    else
        fitMeasure.Dhat(allmonitorsite(monidex)) = fitMeasure.Dhat(allmonitorsite(monidex)) - 0.5*timesize*log(2*pi*tmpSigma2) - 0.5/tmpSigma2*(sum(tmpDelta(1:obj.n1).^2) + ...
            sum((tmpDelta - obj.X0*beta0temp{monidex,1}').^2));
    end
    
    for j0 = 1:(numel( beta1temp{monidex,1}));
        psoutput.thetaCI{monidex,1}{j0} = quantile(Beta1parameter{ComIndex,1}(:,j0),(1-HPDp)/2);
        psoutput.thetaCI{monidex,2}{j0} = quantile(Beta1parameter{ComIndex,1}(:,j0),1-(1-HPDp)/2);
    end
    
    if (numel( beta1temp{monidex,1}) >1)
        tempbeta = Beta1parameter{ComIndex,1};
        for j0 = 1:(numel( beta1temp{monidex,1})- 1)
            tempbetadiff = tempbeta(:,j0+1) - tempbeta(:,j0);
            psoutput.betadifference{monidex,1}{j0} = mean(tempbetadiff);
            [La,Lb] = FindHPDset(tempbetadiff,HPDp,[]);
            psoutput.betadifferenceHPD{monidex,1}{j0} = La;
            psoutput.betadifferenceHPD{monidex,2}{j0} = Lb;
            
            psoutput.betadifferenceCI{monidex,1}{j0} = quantile(tempbetadiff,(1-HPDp)/2);
            psoutput.betadifferenceCI{monidex,2}{j0} = quantile(tempbetadiff,1-(1-HPDp)/2);
            
            if j0 == 1
                cptMeasure.betadiff(allmonitorsite(monidex)) = mean(tempbetadiff);
                cptMeasure.betabwCI(allmonitorsite(monidex)) = mean(psoutput.betadifferenceCI{monidex,2}{j0} - psoutput.betadifferenceCI{monidex,1}{j0});
                cptMeasure.betabwHPD(allmonitorsite(monidex))  = mean(Lb - La);
            end
        end
    end
    
    %     xlabel(['\sigma^2= ', num2str(sigmavaltemp{monidex,1}(1)),',time= ',num2str(TimeInterval(2:end-1)+1)])
    xlabel(['(\sigma_s^2= ', num2str(sigmavaltemp{monidex,1}(1)),')'], 'Fontsize', 16)
    
    %     plot(1:timesize,yvalueupperHPD,'--r'); hold on;
    %     plot(1:timesize,yvaluelowerHPD,'--r');hold on;
    plot(1:timesize,yvalueupperCI,'--r'); hold on;
    plot(1:timesize,yvaluelowerCI,'--r');hold on;
    Rdata.sigma2(monidex) = sigmavaltemp{monidex,1}(1);
    Rdata.yvalueupperCI(:,monidex) = yvalueupperCI';
    Rdata.yvaluelowerCI(:,monidex) = yvaluelowerCI';
    
    %yvalue=[];
    %yvaluec1=[];
    %yvaluec2=[];
    changetime = zeros(size(TimeInterval,2)-1,1);
    tmpyvalue = [];
    tmpsepredT = [];
    
    tmpYs = []; tmpPred = [];
    Rdata.cpt{monidex} = zeros(1,size(TimeInterval,2) - 1);
    for i = 1:size(TimeInterval,2) - 1,
        %timeaxis=[];
        timeaxis = [TimeInterval(i)+1:1:TimeInterval(i+1)];
        changetime(i) = TimeInterval(i+1);
        yvalue = [timeaxis.*beta1temp{monidex,1}(i) + uvaltemp{monidex,1}+ beta0temp{monidex,1}(i)];
        Xmat = [ones(size(timeaxis,2),1) timeaxis'];
        sefit = sqrt(sigmavaltemp{monidex,1}).*sqrt(diag(Xmat*inv(Xmat'*Xmat)*Xmat'))';  % note sigmavaltemp = sigma^2
        sepred = sqrt(sigmavaltemp{monidex,1}).*sqrt(1+diag(Xmat*inv(Xmat'*Xmat)*Xmat'))';
        yvalueprec1 = yvalue-sepred.*tinv(pp/2,size(timeaxis,2)-2);
        yvalueprec2 = yvalue+sepred.*tinv(pp/2,size(timeaxis,2)-2);
        yvaluec1 = yvalue-sefit.*tinv(pp/2,size(timeaxis,2)-2);
        yvaluec2 = yvalue+sefit.*tinv(pp/2,size(timeaxis,2)-2);
        tmpyvalue = [tmpyvalue yvalue];
        tmpsepredT = [tmpsepredT sepred.*tinv(pp/2,size(timeaxis,2)-2)];
        
        if (1<i && i<=size(TimeInterval,2)-1),
            plot([TimeInterval(i)+1, TimeInterval(i)+1],[-90, 90]); hold on;
        else
        end
        Rdata.cpt{monidex}(i) = TimeInterval(i)+1;
        
        %xlabel(['sigma= ', num2str(sigmavaltemp{monidex,1}(1)),',time= ',num2str(TimeInterval(2:end-1))])
        
        %         plot(timeaxis,(yvalue),'r');  % plot the predictive interval
        %         hold on;
        %         plot(timeaxis,(yvalueprec1),'g');
        %         hold on;
        %         plot(timeaxis,(yvalueprec2),'g');
        %         hold on;
        
        tmpYs = [tmpYs, yvalue,yvalueprec1,yvalueprec2];
        tmpPred = [tmpPred, yvalue];
        
        %plot(timeaxis,yvaluec1,'g');
        %hold on;
        %plot(timeaxis,yvaluec2,'g');
        
    end
    tmpYs = [tmpYs, ((1./multipleconst).*alldata(2:end,allmonitorsite(monidex)))', tmpYs2];
    ylim([min(tmpYs), max(tmpYs)+0.06]);
    Rdata.miny(monidex) = min(tmpYs); Rdata.maxy(monidex) = max(tmpYs)+0.06;
    
    tmpPred = tmpPred';
    fitMeasure.Pred(:,allmonitorsite(monidex)) = tmpPred;
    fitMeasure.SSE(allmonitorsite(monidex)) = sum((tmpPred - obsY).^2);
    fitMeasure.MSE(allmonitorsite(monidex)) = sum((tmpPred - obsY).^2)/(timesize - 2*maxtimeIn);
    fitMeasure.BIC(allmonitorsite(monidex)) = log(sum((tmpPred - obsY).^2)/timesize) + 2*maxtimeIn/timesize*log(timesize);
    
    psoutput.yvalue{monidex,1}=tmpyvalue;
    psoutput.sepredT{monidex,1}=tmpsepredT;
    psoutput.changetime{monidex,1}=changetime;
    psoutput.HPDsetdata{monidex,1}=AllHPDsetdata;
end
end


function [Clustergroup,SS,SS2] = makelinkage(xallIf,c)
global sitesize datachoice filename simu Clustergroup0 model rep
modellist = {'SRE (model 1)', 'NSE (model 2)', 'NRE (model 3)'};
if datachoice == 1
    Name=char('Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Washington DC','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky',...
        'Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire	',...
        'New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina',...
        'South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming');
end
Niteration = size(xallIf,1);
PCount = zeros(sitesize,sitesize);
for k = 1:Niteration,
    for j = 1:sitesize,
        jsite = j;
        Nstar = size(xallIf{k,1},1);
        N = [];
        for i = 1:Nstar,
            N(i) = size(xallIf{k,1}{i,1}.Siteincluster,2);% N(i): the number of sites in cluster i.
            if (min(abs(xallIf{k,1}{i,1}.Siteincluster-j))==0)
                for ss=1:size(xallIf{k,1}{i,1}.Siteincluster,2),
                    PCount(j,xallIf{k,1}{i,1}.Siteincluster(ss))=PCount(j,xallIf{k,1}{i,1}.Siteincluster(ss))+1;
                end
            else
            end
        end
    end
end
DissMat= 1-PCount./Niteration;
pdistvec=[];
for i=1:sitesize-1,
    pdistvec = [pdistvec DissMat(i,i+1:end)];
end
clusterRe = linkage(pdistvec,'ward');
%  dendrogram(clusterRe,sitesize, 'LABELS', Name);
if model ~= 3
    [H,T] = dendrogram(clusterRe,sitesize, 'LABELS', Name,'orient','left','colorthreshold',1.5);set(H,'LineWidth',2);axis tight;
else
    [H,T] = dendrogram(clusterRe,sitesize, 'LABELS', Name,'orient','left','colorthreshold',1.3);set(H,'LineWidth',2);axis tight;
end
title(strcat('Model=', modellist(model)),'Fontsize',15)
orient landscape
print('-painters', '-dpsc2', '-r600', '-append', filename)

% get the figure
[H,T] = dendrogram(clusterRe,sitesize, 'LABELS', Name,'orient','left','colorthreshold',1.3);
set(H,'LineWidth',2);axis tight;
title('')
set(gca,'FontSize',13)
orient landscape
print('-painters', '-dpsc2', '-r600', strcat('DendrogramM_',num2str(rep),'.eps'))

ClusterIndex = cluster(clusterRe,'maxclust',c);
Clustergroup = cell(c,1);
for i = 1:sitesize,
    for j = 1:c
        if  ClusterIndex(i)==j,
            Clustergroup{j,1} = [Clustergroup{j,1} i];
        end
    end
end

labs0 = zeros(1,sitesize);
for i = 1:size(Clustergroup);
    labs0(Clustergroup{i}) = i;
end
n1p = 0; n2p = 0;
for i = 1:sitesize
    for j = (i+1):sitesize
        if labs0(i)==labs0(j)
            n1p = n1p+1;
        elseif labs0(i)~=labs0(j)
            n2p = n2p+1;
        end
    end
end

Sensitivity = zeros(Niteration,1); Specificity = zeros(Niteration,1);
for k = 1:Niteration
    labs = zeros(1,sitesize); Nstar = size(xallIf{k,1},1);
    for i = 1:Nstar
        labs(xallIf{k,1}{i,1}.Siteincluster) = i;
    end
    n11 = 0; n22 = 0;
    for i = 1:sitesize
        for j = (i+1):sitesize
            if (labs(i)==labs(j)) && (labs0(i)==labs0(j))
                n11 = n11+1;
            elseif (labs(i)~=labs(j)) && (labs0(i)~=labs0(j))
                n22 = n22+1;
            end
        end
    end
    Sensitivity(k) = n11/n1p; Specificity(k) = n22/n2p;
end
labsC = labs0; 
SS = mean(2-Sensitivity-Specificity); % central
SS2 = []; % true
if simu == 1
    subplot(1,2,1); % iteraton versus central
    scatter(Sensitivity,Specificity,36,'k','filled'); xlim([0,1]); ylim([0,1]); grid on
    xlabel('Sensitivity','FontSize',18); ylabel('Specificity','FontSize',18);
    SensitivityC = Sensitivity; SpecificityC = Specificity;
    
    % iteraton versus true
    labs0 = zeros(1,sitesize);
    for i = 1:size(Clustergroup0);
        labs0(Clustergroup0{i}) = i;
    end
    n1p = 0; n2p = 0;
    for i = 1:sitesize
        for j = (i+1):sitesize
            if labs0(i)==labs0(j)
                n1p = n1p+1;
            elseif labs0(i)~=labs0(j)
                n2p = n2p+1;
            end
        end
    end
    
    Sensitivity = zeros(Niteration,1); Specificity = zeros(Niteration,1);
    for k = 1:Niteration
        labs = zeros(1,sitesize); Nstar = size(xallIf{k,1},1);
        for i = 1:Nstar
            labs(xallIf{k,1}{i,1}.Siteincluster) = i;
        end
        n11 = 0; n22 = 0;
        for i = 1:sitesize
            for j = (i+1):sitesize
                if (labs(i)==labs(j)) && (labs0(i)==labs0(j))
                    n11 = n11+1;
                elseif (labs(i)~=labs(j)) && (labs0(i)~=labs0(j))
                    n22 = n22+1;
                end
            end
        end
        Sensitivity(k) = n11/n1p; Specificity(k) = n22/n2p;
    end
    SS2 = mean(2-Sensitivity-Specificity);
    subplot(1,2,2);
    scatter(Sensitivity,Specificity,36,'k','filled'); xlim([0,1]); ylim([0,1]); grid on
    hold on;
    % get SS for C versus T
    n11 = 0; n22 = 0; 
    for i = 1:sitesize
        for j = (i+1):sitesize
            if (labsC(i)==labsC(j)) && (labs0(i)==labs0(j))
                n11 = n11+1;
            elseif (labsC(i)~=labsC(j)) && (labs0(i)~=labs0(j))
                n22 = n22+1;
            end
        end
    end
    SensitivityCT = n11/n1p; SpecificityCT = n22/n2p; 
    labsT = labs0;
else
    scatter(Sensitivity,Specificity,36,'ko'); xlim([0,1]); ylim([0,1]); grid on % central
    %title(strcat('Model=', modellist(model)),'Fontsize',15)
    SS2 = []; %mean(2-Sensitivity-Specificity);
end
if simu == 0
    save('SS.mat','Sensitivity','Specificity') % central
else
    save(strcat('SS_',num2str(rep),'.mat'),'SensitivityC','SpecificityC','Sensitivity',...
        'Specificity','labsC','labsT','SensitivityCT','SpecificityCT') % central, true
end
xlabel('Sensitivity','FontSize',18); ylabel('Specificity','FontSize',18);
title(strcat('max d = ',num2str(c)),'FontSize',18)
set(gca,'FontSize',16)
orient landscape
print('-painters', '-dpsc2', '-r600', 'SSc_real.eps')
end

function [out] = getDistance2(data, Clustergroup0, Clustergroup1)
global dataPath rep model MinNtimesw timesize sitesize simu 
N = floor(timesize/MinNtimesw);
Name=char('Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Washington DC','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky',...
    'Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire	',...
    'New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina',...
    'South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming');

computP = 1;
if computP == 1
    if simu == 0
        % load(strcat(dataPath,'Ressandbeta',num2str(model),'_',num2str(rep),'.mat'))
        load('Ressandbeta.mat')
    else % for simulation I am going to calculate it here
        Ressval = cell(N,sitesize);
        beta0s = cell(N,sitesize);
        beta1s = cell(N,sitesize);
        for j = 1:sitesize
            for i=1:N % the number of time intervals.
                if i==1 %only one time inveral in the site j
                    sitedata = data(:,j); times = 1:timesize;
                    [Ressval{i,j},beta1s{i,j},beta0s{i,j}] = Ress(sitedata',times);
                else
                    sitedata = data(:,j);
                    B = Combination(i,timesize-i*MinNtimesw);
                    m = size(B,1); k = size(B,2);
                    times = 1:timesize;
                    celly = cell(m,k);
                    timesinterval=cell(m,k);
                    for l=1:m
                        timesinterval(l,:) = mat2cell(times,1,B(l,:));
                        celly(l,:)= mat2cell(sitedata',1,B(l,:));
                    end
                    [Ressval{i,j},beta1s{i,j},beta0s{i,j}] = cellfun(@Ress,celly,timesinterval);
                end
            end
        end
    end
    Dist = zeros(sitesize); count = 1; times = 1:timesize;
    Info = cell(sitesize);
    Hcps = cell(1,sitesize); Hbeta0s = cell(1,sitesize); Hbeta1s = cell(1,sitesize); HRa = zeros(1,sitesize);
    for s1 = 1:sitesize % find changepoint measures
        D0 = 1e4;
        for i = 1:N
            if i == 1
                D = timesize*log(Ressval{i,s1}) + (2*i+1)*log(timesize);
                if D<D0
                    D0 = D; Hcps{s1} = timesize; Hbeta0s{s1} = beta0s{i,s1}; Hbeta1s{s1} = beta1s{i,s1};
                end
            else
                B = Combination(i,timesize-i*MinNtimesw);
                m = size(B,1); k = size(B,2); D = 0;
                for l = 1:m
                    D = timesize*log(sum(Ressval{i,s1}(l,:))) + (2*i+1)*log(timesize);
                    if D<D0
                        D0 = D; Hcps{s1} = B(l,:); Hbeta0s{s1} = beta0s{i,s1}(l,:); Hbeta1s{s1} = beta1s{i,s1}(l,:);
                    end
                end
            end
        end
        if simu == 1 % detect accuracy
            load('true_cpt.mat')
            for r = 1:length(Clustergroup0)
                if any(Clustergroup0{r} == s1)
                  HRa(s1) = 0;
                  if length(Hcps{s1}) > 1
                      tmp0 = cumsum([1 Hcps{s1}]);
                      obs_cpt = tmp0(2:(length(tmp0)-1)); lens = length(obs_cpt);
                      for k = 1:lens
                          if any(true_cpt{r} == obs_cpt(k))
                              HRa(s1) = HRa(s1)+1;
                          end
                      end
                      if lens > length(true_cpt{r})
                          HRa(s1) = HRa(s1) - (lens-length(true_cpt{r}));
                      end
                      HRa(s1) = HRa(s1)/length(true_cpt{r});
                  end
                end
            end
        end
    end
    out.HRa = HRa; out.Hbeta0s = Hbeta0s; out.Hbeta1s = Hbeta1s; out.Hcps = Hcps;
    for s1 = 2:sitesize % find clustering measures
        for s2 = 1:(s1-1)
            fprintf('%6d', count); count = count+1;
            if(~mod(count,20))
                fprintf('\n')
            end
            D0 = 1e4; sites = [s1,s2];
            for i=1:N % the number of time intervals.
                if i==1 %only one time inveral in the site j
                    D = 0;
                    betamle = mean([beta1s{i,s1}, beta1s{i,s2}]);
                    for j = 1:2
                        sitedata = data(:,sites(j));
                        D = D + timesize*log(sum(( sitedata-mean(sitedata) - betamle*(times'-mean(times')) ).^2)/Ressval{i,sites(j)});
                    end
                    D = D-i*log(timesize);
                    if D<D0
                        D0 = D; Info{s1,s2} = [i-1,0,betamle];
                    end
                else
                    B = Combination(i,timesize-i*MinNtimesw);
                    m = size(B,1); k = size(B,2); D = 0;
                    for l=1:m
                        betamle = zeros(1,k);
                        for p = 1:k
                            betamle(p) = mean([beta1s{i,s1}(l,p), beta1s{i,s2}(l,p)]);
                        end
                        timesinterval = mat2cell(times,1,B(l,:));
                        for j = 1:2
                            sitedata = data(:,sites(j));
                            celly = mat2cell(sitedata',1,B(l,:));
                            res = 0;
                            for p = 1:k
                                res = res+sum(( celly{p}-mean(celly{p}) - betamle(p)*(timesinterval{p}-mean(timesinterval{p})) ).^2);
                            end
                            D = D + timesize*log(res/sum(Ressval{i,sites(j)}(l,:)));
                        end
                        D = D-i*log(timesize);
                        if D<D0
                            D0 = D; Info{s1,s2} = [i-1,B(l,:),betamle];
                        end
                    end
                end
            end
            Dist(s1,s2) = D0;
        end
    end
    for s1 = 2:sitesize
        for s2 = 1:(s1-1)
            Dist(s2,s1) = Dist(s1,s2);
        end
    end
    pdistvec=[];
    for i=1:sitesize-1,
        pdistvec = [pdistvec Dist(i,i+1:end)];
    end
    Dmin = min(pdistvec); Dmax = max(pdistvec);
    pdistvec = (pdistvec-Dmin)/(Dmax-Dmin);
    clusterRe = linkage(pdistvec, 'ward');
    %  dendrogram(clusterRe,sitesize, 'LABELS', Name);
    [H,T] = dendrogram(clusterRe,sitesize, 'LABELS', Name,'orient','left','colorthreshold',1.3);
    set(H,'LineWidth',2);axis tight;
    title('')
    set(gca,'FontSize',13)
    orient landscape
    print('-painters', '-dpsc2', '-r600', strcat('DendrogramP_',num2str(rep),'.eps'))
    out.Dist = Dist;
    out.Info = Info;
end

% for true labeling
if simu == 1
    cvec = 3:8; clen = length(cvec);
    Sensitivity = zeros(clen,1); Specificity = zeros(clen,1);
    for cind = 1:clen
        c = cvec(cind);
        ClusterIndex = cluster(clusterRe,'maxclust',c);
        Clustergroup = cell(c,1);
        for i = 1:sitesize,
            for j = 1:c
                if  ClusterIndex(i)==j,
                    Clustergroup{j,1} = [Clustergroup{j,1} i];
                end
            end
        end
        
        labs0 = zeros(1,sitesize);
        for i = 1:size(Clustergroup0);
            labs0(Clustergroup0{i}) = i;
        end
        n1p = 0; n2p = 0;
        for i = 1:sitesize
            for j = (i+1):sitesize
                if labs0(i)==labs0(j)
                    n1p = n1p+1;
                elseif labs0(i)~=labs0(j)
                    n2p = n2p+1;
                end
            end
        end
        labs = zeros(1,sitesize);
        for i = 1:size(Clustergroup);
            labs(Clustergroup{i}) = i;
        end
        n11 = 0; n22 = 0;
        for i = 1:sitesize
            for j = (i+1):sitesize
                if (labs(i)==labs(j)) && (labs0(i)==labs0(j))
                    n11 = n11+1;
                elseif (labs(i)~=labs(j)) && (labs0(i)~=labs0(j))
                    n22 = n22+1;
                end
            end
        end
        Sensitivity(cind) = n11/n1p; Specificity(cind) = n22/n2p;
        text(Sensitivity(cind),Specificity(cind),num2str(c),...
            'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7],'FontSize',11)
    end
    out.SS2 = [cvec', Sensitivity, Specificity, 2-Sensitivity-Specificity];
end

% for central labeling
cvec = 3:8; clen = length(cvec);
Sensitivity = zeros(clen,1); Specificity = zeros(clen,1);
for cind = 1:clen
    c = cvec(cind);
    ClusterIndex = cluster(clusterRe,'maxclust',c);
    Clustergroup = cell(c,1);
    for i = 1:sitesize,
        for j = 1:c
            if  ClusterIndex(i)==j,
                Clustergroup{j,1} = [Clustergroup{j,1} i];
            end
        end
    end
    
    labs0 = zeros(1,sitesize);
    for i = 1:size(Clustergroup1);
        labs0(Clustergroup1{i}) = i;
    end
    n1p = 0; n2p = 0;
    for i = 1:sitesize
        for j = (i+1):sitesize
            if labs0(i)==labs0(j)
                n1p = n1p+1;
            elseif labs0(i)~=labs0(j)
                n2p = n2p+1;
            end
        end
    end
    labs = zeros(1,sitesize);
    for i = 1:size(Clustergroup);
        labs(Clustergroup{i}) = i;
    end
    n11 = 0; n22 = 0;
    for i = 1:sitesize
        for j = (i+1):sitesize
            if (labs(i)==labs(j)) && (labs0(i)==labs0(j))
                n11 = n11+1;
            elseif (labs(i)~=labs(j)) && (labs0(i)~=labs0(j))
                n22 = n22+1;
            end
        end
    end
    Sensitivity(cind) = n11/n1p; Specificity(cind) = n22/n2p;
    if simu~=1
        text(Sensitivity(cind),Specificity(cind),num2str(c),...
            'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7],'FontSize',11)
    else
        subplot(1,2,1); hold on; text(Sensitivity(cind),Specificity(cind),num2str(c),...
            'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7],'FontSize',11)
    end
end
out.SS = [cvec', Sensitivity, Specificity, 2-Sensitivity-Specificity];
end

function [out] = computefromTable(Clustergroup, Clustergroup0)
global sitesize

len = sitesize * (sitesize - 1) / 2;
rul = zeros(len,2); k = 1;
for i = 1:sitesize
    for j = (i+1):sitesize
        rul(k,:) = [i, j]; k = k+1;
    end
end
lab0 = zeros(1,len); lab1 = zeros(1,len);

for i = 1:size(Clustergroup0,1)
    for i1 = 1:length(Clustergroup0{i})
        for i2 = (i1+1):length(Clustergroup0{i})
            s1 = min(Clustergroup0{i}([i1 i2]));
            s2 = max(Clustergroup0{i}([i1 i2]));
            for k = 1:len
                if rul(k,1) == s1 && rul(k,2) == s2
                    lab0(k) = 1;
                end
            end
        end
    end
end

for i = 1:size(Clustergroup,1)
    for i1 = 1:length(Clustergroup{i})
        for i2 = (i1+1):length(Clustergroup{i})
            s1 = min(Clustergroup{i}([i1 i2]));
            s2 = max(Clustergroup{i}([i1 i2]));
            for k = 1:len
                if rul(k,1) == s1 && rul(k,2) == s2
                    lab1(k) = 1;
                end
            end
        end
    end
end

tab = crosstab(lab0, lab1); n = sum(sum(tab));
n11 = tab(1,1); n22 = tab(2,2); n12 = tab(1,2); n21 = tab(2,1);
n1p = sum(tab(1,:)); n2p = sum(tab(2,:));
np1 = sum(tab(:,1)); np2 = sum(tab(:,2));

out.tab = tab;
out.chisq = (n11 - n1p)^2/n1p + (n22 - n2p)^2/n2p;
out.phi = 1 - out.chisq/n;
out.yuleQ = (n11*n22 - n12*n21) / (n11*n22 + n12*n21);
out.specificity = n11/n1p;   % 1,1 = both are 0
out.sensitivity = n22/n2p;
pa = (n11+n22)/n; pe = np1/n*n1p/n + np2/n*n2p/n;
out.Cohenkappa = (pa - pe) / (1 - pe);
end

function [LBout,UBout] = FindHPDset(Samples,p,npoints)
%Function to find the 100p % HPD set based on Samples
if isempty(npoints)
    npoints=200;
end

[f,x] = ksdensity(Samples,'npoints',npoints); N = size(Samples,1); maxf = max(f);
step = maxf/npoints;
HPDdensity = (step:step:maxf);
NHPD = size(HPDdensity,2);
LB = cell(1,NHPD); UB = cell(1,NHPD); Prob = zeros(1,NHPD);
for i=1:NHPD
    indices0 = find(HPDdensity(NHPD-i+1) < f);
    if ~isempty(indices0)
        indices1 = find(diff(indices0)> 1);
        if isempty(indices1)
            LB{i} = x(indices0(1)); UB{i} = x(indices0(end));
        elseif (size(indices1,1)==1)
            LB{i} = [x(indices0(1)) x(indices0(indices1(1)+1))];
            UB{i} = [x(indices0(indices1(1))) x(indices0(end))];
        else
            LB{i} = x(indices0(1)); UB{i} = [];
            for j=1:(size(indices1,2)-1)
                LB{i} = [LB{i} x(indices0(indices1(j)+1))];
                UB{i} = [UB{i} x(indices0(indices1(j)))];
            end
            UB{i} =[UB{i} x(indices0(end))];
        end
    end
    Ns = size(LB{i},2);
    count = zeros(1,Ns);
    for j=1:Ns
        count(j) = sum((LB{i}(j) <= Samples).*(Samples <= UB{i}(j)));
    end
    Prob(i) = sum(count)/N;
end
[minval indexmin] = min(abs(Prob - p));
LBout = LB{indexmin};
UBout = UB{indexmin};
end

function [R,neff,V,W,B] = psrf(varargin)
%PSRF Potential Scale Reduction Factor
%
%   [R,NEFF,V,W,B] = PSRF(X) or
%   [R,NEFF,V,W,B] = PSRF(x1,x2,...,xs)
%   returns "Potential Scale Reduction Factor" (PSRF) for
%   collection of MCMC-simulations. X is a NxDxM matrix
%   which contains M MCMC simulations of length N, each with
%   dimension D. MCMC-simulations can be given as separate
%   arguments x1,x2,... which should have the same length.
%
%   Returns
%     R     PSRF in a row vector of length D
%     neff  estimated effective number of samples M*N*V/B
%     V     estimated mixture-of-sequences variances
%     W     estimated within sequence variances
%     B     estimated between sequence variances
%
%   The idea of the PSRF is that if R is not near 1 (below 1.1 for
%   example) one may conclude that the tested samples were not from
%   the same distribution (chain might not have been converged
%   yet).
%
%   If only one simulation is given, the factor is calculated
%   between first and last third of the chain. Note that use of
%   only one chain will produce over-optimistic result.
%
%   Method is from:
%      Brooks, S.P. and Gelman, A. (1998) General methods for
%      monitoring convergence of iterative simulations. Journal of
%      Computational and Graphical Statistics. 7, 434-455. Note that
%      this function returns square-root definiton of R (see Gelman
%      et al (2003), Bayesian Data Analsyis p. 297).
%
%   See also
%     CPSRF, MPSRF, IPSRF

% Copyright (C) 1999 Simo Särkk?E
% Copyright (C) 2003 Aki Vehtari
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

% 2004-01-22 Aki.Vehtari@hut.fi Added neff, R^2->R, and cleaning

% In case of one argument split to two halves (first and last thirds)
onechain=0;
if nargin==1
    X = varargin{1};
    if size(X,3)==1
        n = floor(size(X,1)/3);
        x = zeros([n size(X,2) 2]);
        x(:,:,1) = X(1:n,:);
        x(:,:,2) = X((end-n+1):end,:);
        X = x;
        onechain=1;
    end
elseif nargin==0
    error('Cannot calculate PSRF of scalar');
else
    X = zeros([size(varargin{1}) nargin]);
    for i=1:nargin
        X(:,:,i) = varargin{i};
    end
end

[N,D,M]=size(X);

if N<1
    error('Too few samples');
end

% Calculate means W of the variances
W = zeros(1,D);
for n=1:M
    x = X(:,:,n) - repmat(mean(X(:,:,n)),N,1);
    W = W + sum(x.*x);
end
W = W / ((N-1) * M);

% Calculate variances B (in fact B/n) of the means.
Bpn = zeros(1,D);
m = mean(reshape(mean(X),D,M)');
for n=1:M
    x = mean(X(:,:,n)) - m;
    Bpn = Bpn + x.*x;
end
Bpn = Bpn / (M-1);

% Calculate reduction factors
S = (N-1)/N * W + Bpn;
R = (M+1)/M * S ./ W - (N-1)/M/N;
V = R .* W;
R = sqrt(R);
B = Bpn*N;
neff = min(M*N*V./B,M*N);
if onechain & (nargout>1)
    neff=neff*3/2;
end
end

function [B] = Combination(Kstar,S)
global MinNtimesw
global timesize
if Kstar==1,
    B=timesize;
else
    args = cell(1,Kstar);
    for i=1:Kstar
        args{i} = [0:1:S];
    end
    [A{1:Kstar}] = ndgrid(args{1:Kstar});
    A = reshape(cat(Kstar+1,A{:}),(S+1)^Kstar,Kstar);
    C1 = (sum(A,2) == S);
    B = A(C1,:)+MinNtimesw ;
end
end

function [obj] = computefromInt(NpointinTimeIn)
reducedX0 = 1;
global timesize inc_ints
p = length(NpointinTimeIn); % number of intervals = Ks + 1
if p == 1
    reducedX0 = 0;
end
Tvec = cumsum([1 NpointinTimeIn]); Tvec(end) = timesize;
obj.X0 = zeros(timesize, p); obj.X1 = zeros(timesize, p);
for i = 1:p
    vec = Tvec(i):(Tvec(i+1)-1);
    obj.X0(vec,i) = 1;
    obj.X1(vec,i) = vec;
end
obj.X0(timesize, p) = 1;
% if(reducedX0)
%     obj.X0 = obj.X0(Tvec(2):end,2:end);
% else
%     obj.X0 = 1;
% end
obj.X1(timesize, p) = timesize;
obj.D = obj.X0'*obj.X0;
obj.n1 = NpointinTimeIn(1);
if inc_ints ~= 1
    obj.Sigma = eye(timesize);
    if(reducedX0)
        obj.Sigma((obj.n1+1):end,(obj.n1+1):end) = eye(size(obj.X0,1)) - obj.X0*inv(obj.D)*obj.X0';
    end
else
    obj.Sigma = eye(size(obj.X0,1)) - obj.X0*inv(obj.D)*obj.X0';
end
end

function [RSS,beta1,beta0]=Ress(datavalue,timeIn)
x=[ones(size(timeIn,2),1) timeIn'];
[betas,bint,Ressval] =regress(datavalue',x);
beta1=betas(2);
beta0=betas(1);
RSS=sum(Ressval.^2);
end