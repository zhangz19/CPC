function [] = cpc(ID)
% by Zhen Zhang (zhangquake@gmail.com), 09/22/2013
global W M timesize sitesize MinNtimesw kstar Nint multipleconst
global SumHpart1addLogMulnomial Hpart2 Hinfo %#ok<*NUSED>
global alphalambda betalambda alphasig betasig simu alphatau2 betatau2
global priorchoice delta2 datachoice modelchoice random
global integrated useT inc_int inc_ints dataPath aAlpha0 bAlpha0 fixsmalllambda
global loglike0 gap gammas

% random = 1; %0: without random effect
% modelchoice = 1; %1: spatial, 2: hetero, 3: random effect

simu = 1; %run simulation, turn it = 0 for real data
verbose = 0;
datachoice = 1; %1: U.S. state, %2: iowa state, %3: Yeast cell data
% totIter = 15e3; burn = 1e4;
multipleconst = 1;
integrated = 1;
inc_int = 0; % inc_int=1: include intercept alpha_r 
inc_ints = 1; % inc_int=1: include intercept alpha_s
updtAlpha0 = 1;
fixsmalllambda = 0;
dataPath = './';

ch = str2double(num2str(ID));   % if simu == 1: ch = 1- 30*6 = 180
if simu == 0
    totIter = 6e3; burn = 5e3;
    nChain = 4; % 2 if simu == 1
    useT = 0;
else
    totIter = 6e3; burn = 5e3;
    %totIter = 3e3; burn = 2e3;
    nChain = 6;
    useT = 0;
end
ch0 = ch;
rep = ceil(ch/nChain);
ch = ch - (rep-1)*nChain;
inits = ch;
ch = 3; % ch now represents the id of model

% rng('default'); rng(ch0*210); % set random seed
rng('default'); rng(ch0*24);

switch ch
    case 1
        random = 1; modelchoice = 1;
    case 2
        random = 1; modelchoice = 3;
    case 3
        random = 0; modelchoice = 3;
end

priorchoice = 1; % 1: flat, 2: informative diffuse prior for the intercept
delta2 = 0.8^2;  % p(alpha) ~ N(0,delta2*I)

% meanlambda = 1e-3; varlambda = 1e-270;
meanlambda = 2; varlambda = 1/.001;
betalambda = varlambda/meanlambda;
alphalambda = meanlambda/betalambda; %for lambda

meantau2 = 0.01; vartau2 = 10^2;
alphatau2 = 2+meantau2^2/vartau2;
betatau2 = 1/(meantau2*(alphatau2-1)); %[1/((alphasig-1)*betasig) 1/((alphasig-1)^2*(alphasig-2)*betasig^2)]
% alphatau2 = 2; betatau2 = 1;

meansigma2 = 0.01; varsigma2 = 10^2;
alphasig = 2+meansigma2^2/varsigma2;
betasig = 1/(meansigma2*(alphasig-1));
% alphasig = 2; betasig = 1;

alpha = 1e2;
aAlpha0 = 2; bAlpha0 =  1/.001;

if datachoice == 1
    load('neighbor.txt')
else
    load('iowa_neighbor.mat'); neighbor = NB; %larger
end

if simu == 1
    %     SimulateData(datachoice);
    load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat'))
    MinNtimesw = 7;
    %     GenerateRSSandbetamodel1(data, modelchoice); %#ok<NODEF>
    data = data(2:end,:);
else
    if datachoice == 1
        load('data.txt') % for real data
        MinNtimesw = 7;
        %     GenerateRSSandbetamodel1(data, modelchoice); %#ok<NODEF>
        data = data(2:end,:);
    elseif  datachoice == 2
        data = load('iowa09-12.txt'); data = data'; data = data(5:35,:);
        MinNtimesw = 7;
    elseif  datachoice == 3
        data = load('yeast cell data.txt');
        MinNtimesw = 3;
    end
end

[timesize,sitesize] = size(data);
W = neighbor;
M = diag(sum(W,1));

if datachoice == 3
    M = eye(sitesize); W = eye(sitesize);  % non-spatial
end

eigs = eig(sqrt(inv(M))*W*sqrt(inv(M)));
lgamma = max(1/min(eigs),-1); ugamma = 1/max(eigs);
gap = 1e-3;
gammas = (lgamma+gap):gap:(ugamma-gap); len = length(gammas);
loglike0 = zeros(1,len);
for i = 1:len
    gamma = gammas(i);
    B = (M-gamma*W);
    loglike0(i) = 0.5*sum(log(eig(B)));
end

kstar = floor((timesize-1)/MinNtimesw-1);
Nint = kstar + 1; %floor(timesize/MinNtimesw);

Createdata(data, ch, inits, rep);
Precalculation(ch, inits, rep);

load(strcat('InitialValue',num2str(ch),'_',num2str(inits),'_',num2str(rep),'.mat'))
load(strcat('outPrecalculation',num2str(ch),'_',num2str(inits),'_',num2str(rep),'.mat'))

bloop = 12;
% rand('state',bloop*ch0+rep*30*nChain);
% randn('state',bloop*ch0+rep*30*nChain+2^20);

x{1}.alpha = alpha;

xall = cell(1,totIter-burn);
tic
for ss = 1:totIter
    if verbose == 1
        display(ss)
        tmpN = [];tmpNS = [];
        for NC = 1:size(x,1)
            tmpN = [tmpN x{NC,1}.NTimeIntervals];
            tmpNS = [tmpNS numel(x{NC,1}.Siteincluster)];
        end
        display([[tmpN x{1}.tau2 x{1}.lambda];[tmpNS x{1}.gamma x{1}.alpha]])
    else
        fprintf('%6d', ss)
        if ~mod(ss,20)
            fprintf('\n')
        end
    end
    
    if ss > burn   % record after burning time
        xall{ss-burn} = x;
    end
    
    [x] = UpdatesiteNofixedsigma(x);  %#ok<*NASGU>
    [x] = Updatesigma2Nofixedsigma(x);
    [x] = UpdateClusterInfoNofixedsigma(x);
    [x] = UpdatelambdaNofixedsigma(x);
    if random == 1
        [x] = UpdateU(x);
        if modelchoice == 1
            [x] = UpdateGamma(x);
        end
        [x] = UpdateTau2(x);
    end
    if updtAlpha0 == 1
        [x] = UpdateAlpha0(x);
    end
end
CPUtime = toc; CPUtime = CPUtime/60;
fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', totIter, CPUtime)
save(strcat('MCMCoutputCh',num2str(ch),'_ini',num2str(inits),'_',num2str(rep),'.mat'),'xall','CPUtime')

end

function [] = Createdata(data, ch, inits, rep)
global multipleconst timesize sitesize M W simu datachoice
global random modelchoice useT u0 MinNtimesw Lungdata inc_ints inc_int

lambda = 1e-3;
gamma = (inits-2)/3;
tau2 = 0.001;   % empirical

if modelchoice ~= 1
    gamma = 0.0;
end

invM = inv(M);
DU = tau2*inv(eye(size(invM, 1)) - gamma*invM*W)*invM;
if datachoice ~= 3
    randn('state', 200*inits); Uvec = mvnrnd(zeros(1,sitesize),DU);
else
    Uvec = zeros(1,sitesize);
end

if useT == 1
    data0 = data;
    load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat'))
    data = data0;
    Uvec = true_uvec + (inits-1)/1e4;
    tau2 = true_tau2 + (inits-1)/1e4;
    gamma = true_gamma - (inits-1)*0.3;
end

if simu == 1
    x{1,1}.NpointinTimeIn = [10 7 7 7 7];
    x{2,1}.NpointinTimeIn = [7 13 18];
    x{3,1}.NpointinTimeIn = [7 10 11 10];
    x{4,1}.NpointinTimeIn = [17 21];
    x{5,1}.NpointinTimeIn = [11 15 12];
    x{6,1}.NpointinTimeIn = [20 18];
    
    for r = 1:6 % modified 9-22-2013
        x{r,:}.NpointinTimeIn = timesize;
    end
else
    if datachoice == 1
        x{1,1}.NpointinTimeIn = [18 20];
        x{2,1}.NpointinTimeIn = [19 19];
        x{3,1}.NpointinTimeIn = [21 17];
        x{4,1}.NpointinTimeIn = [11 11 16];
    else
        x{1,1}.NpointinTimeIn = timesize;
    end
end

if simu == 1
    if datachoice == 1
        inx = randperm(sitesize); % modified 9-22-2013, need varying seeds
        x{1,:}.Siteincluster= inx(1:9);
        x{2,:}.Siteincluster= inx(10:16);
        x{3,:}.Siteincluster= inx(17:23);
        x{4,:}.Siteincluster= inx(24:34);
        x{5,:}.Siteincluster= inx(35:41);
        x{6,:}.Siteincluster= inx(42:sitesize);
        
    elseif datachoice == 2
        x{1,:}.Siteincluster=1:30;
        x{2,:}.Siteincluster=31:50;
        x{3,:}.Siteincluster=51:70;
        x{4,:}.Siteincluster=71:sitesize;
    end
else
    if datachoice == 1
        % results based on previous running for real data
        x{1,:}.Siteincluster = [1   ,  3   , 16  ,  23  ,  32 ,   33   , 40  ,  41];
        x{2,:}.Siteincluster = [11  ,  13  ,  14   , 15 ,   18  ,  22  ,  24  ,  25 ,   26  ,  35 ,   39   , 47 ,   48];
        x{3,:}.Siteincluster = [2  ,   4  ,   8  ,   9  ,  19 ,   28  ,  29   , 30  ,  31  ,  38  ,  42  ,  43];
        x{4,:}.Siteincluster = [5   ,  6  ,   7  ,  10 ,   12  ,  17  ,  20   , 21 ,   27  ,  34  ,  36  ,  37   , 44   , 45  ,  46 ,   49];
    else
        x{1,:}.Siteincluster=1:sitesize;
    end
end

if(simu)
    Lungdata = data;
else
    if datachoice == 1
        Lungdata = multipleconst*log(data); %log rate
    elseif datachoice == 2
        Lungdata = multipleconst*data;
    elseif datachoice == 3
        Lungdata = multipleconst*data;
    end
end

% now determine the prior mean of CAR model
u0 = zeros(1,sitesize);
if simu ~= 1 && inc_ints ~= 1  % for real data
    for j = 1:sitesize
        ay = Lungdata(1:MinNtimesw,j); ax = 1:MinNtimesw;
        betas =regress(ay, [ones(length(ax),1), ax']);
        u0(j) = betas(1);   % mean vector of random effects
    end
    Lungdata = Lungdata - repmat(u0, [timesize 1]);
    u0 = zeros(1,sitesize);
end

if random == 0
    Uvec = u0;
end

% plot(Lungdata)
ClusterN = size(x,1);

for i = 1:ClusterN,
    datatemp = [];
    xR = [];
    for j = 1:size(x{i,:}.Siteincluster,2),
        xtemp=[];
        datatemp = [datatemp Lungdata(:,x{i,:}.Siteincluster(j))];
        xtemp = [ones(timesize,1) [1:1:timesize]'];
        xR = [xR xtemp];
    end
    yvalue = mat2cell(datatemp,x{i,1}.NpointinTimeIn,size(x{i,:}.Siteincluster,2));
    XRs = mat2cell(xR,x{i,1}.NpointinTimeIn,size(x{i,:}.Siteincluster,2)*2);
    yallvalues = cell(size(yvalue,1),1);
    Xmats = cell(size(yvalue,1),1); labels = cell(1,size(yvalue,1));
    for s = 1:size(yvalue,1),
        yvalues = [];
        XXRs = []; labs = [];
        for j = 1:size(yvalue{s,1},2),
            yvalues = [yvalues;yvalue{s,1}(:,j)];
            XXRs = [XXRs;XRs{s,1}(:,2*j-1:2*j)];
            labs = [labs, j*ones(1, length(yvalue{s,1}(:,j)))];
        end
        yallvalues{s,:} = yvalues;
        Xmats{s,:} = XXRs;  labels{s} = labs;
    end
    out = cellfun(@RessCreatedata,yallvalues,Xmats)'; %+rand(1,size(yvalue,1));
    TMPvec = zeros(1,length(out)); TMPvec0 = zeros(1,length(out));
    for i2 = 1:length(out)
        TMPvec(i2) = out(i2).betaone;  TMPvec0(i2) = out(i2).betazero;
    end
    x{i,1}.parameter = TMPvec; % + rand(1, length(TMPvec))*abs(mean(TMPvec))/5000*inits;  % MAY BE CHANGED
    if inc_int == 1
        x{i,1}.beta0 = TMPvec0; % + rand(1, length(TMPvec0))*abs(mean(TMPvec0))/5000*inits;
    elseif inc_ints == 1
        for j = 1:1:size(x{i,:}.Siteincluster,2)
            x{i,1}.beta0{j} = TMPvec0;
        end
    end
    
    if useT == 1
        x{i,1}.parameter = true_bvec{i}'+ rand(1, length(true_bvec{i}))*abs(mean(true_bvec{i}))/100*inits;
        if inc_ints ~= 1
            x{i,1}.beta0 = true_avec{i}';
        else
            for j = 1:1:size(x{i,:}.Siteincluster,2)
                x{i,1}.beta0{j} = true_avec{i}';
            end
        end
    end
    
    for j = 1:size(x{i,:}.Siteincluster,2)
        tmpres = [];
        for i2 = 1:length(out)
            tmpres = [tmpres, out(i2).Ressval(labels{i2} == j)'];
        end
        x{i,1}.sigma2(j) = var(tmpres);
        if useT == 1
            x{i,1}.sigma2(j) = true_sigma2(x{i,:}.Siteincluster(j));
        end
    end
end

for k = 1:ClusterN,
    datatemp = [];
    data = []; Us = [];
    sn=size(x{k,:}.Siteincluster,2);%the number of the site in the cluster
    x{k,1}.NTimeIntervals = size(x{k,1}.NpointinTimeIn,2);
    for j=1:sn,
        tmpU = Uvec(x{k,:}.Siteincluster(j));
        Us = [Us tmpU];
        %datatemp = [datatemp Lungdata(:,x{k,:}.Siteincluster(j)) - tmpU];
        datatemp = [datatemp Lungdata(:,x{k,:}.Siteincluster(j))];
    end
    stdtemp = std(datatemp);
    data = mat2cell(datatemp,x{k,1}.NpointinTimeIn,sn);
    if x{k,:}.NTimeIntervals == 1,
        x{k,1}.NpointinTimeIn(1) = timesize;
        for i = 1:sn,
            for s = 1:x{k,:}.NTimeIntervals,
                x{k,1}.data{i,1}.TimeIntervals{s,:} = data{s,1}(:,i)';
            end
            %             x{k,1}.sigma2(i) = stdtemp(i)*stdtemp(i);
        end
    else
        for i = 1:sn,
            for s = 1:x{k,:}.NTimeIntervals,
                x{k,1}.data{i,1}.TimeIntervals{s,:} = data{s,1}(:,i)';
            end
            %             x{k,1}.sigma2(i) = stdtemp(i)*stdtemp(i);
        end
    end
    x{k,1}.u = Us;
    x{k,1}.lambda = lambda;
    x{k,1}.tau2 = tau2;
    x{k,1}.gamma = gamma;
end

filename = strcat('InitialValue',num2str(ch),'_',num2str(inits),'_',num2str(rep),'.mat');

save(filename,'x');
end

function [] = Precalculation(ch, inits, rep)
global MinNtimesw timesize Nint integrated alphasig betasig inc_int

SumHpart1addLogMulnomial = cell(1,Nint);
Hpart2 = cell(1,Nint);
Hinfo = cell(1,Nint);

for i = 1:Nint % the number of time intervals.
        B = Combination(i,timesize-i*MinNtimesw);
        m = size(B,1); k = size(B,2);
        Hvalue1 = zeros(m,1);
        for l = 1:m
            Hvalue1(l,1) = getHpart1(B(l,:),integrated,inc_int,alphasig,betasig);
            Hpart2{i}{l} = getHpart2(B(l,:), integrated,inc_int,betasig);
            Hinfo{i}{l} = getHpartInfo(B(l,:));
        end
        SumHpart1addLogMulnomial{i} = Hvalue1 + logmultinomial(B-MinNtimesw,ones(1,k)./i);
end
save(strcat('outPrecalculation',num2str(ch),'_',num2str(inits),'_',num2str(rep),'.mat'),'SumHpart1addLogMulnomial','Hpart2','Hinfo');
end

function [Hvalpart1] = getHpart1(NpointinTimeIn,integrated,inc_int,alphasig,betasig)
global timesize
obj = computefromInt(NpointinTimeIn);
p = length(NpointinTimeIn); % number of intervals = Ks + 1
fac = 0.5*(timesize - 2*(p-1)-2);
Hvalpart1 = - fac*log(2*pi) - 0.5*log(det(obj.X0'*obj.X0))...
    - 0.5*log(det(obj.X1'*obj.Sigma*obj.X1));
if integrated == 1
    nstar = 2*alphasig + timesize - 2*(p-1) - 2;
    Hvalpart1 = Hvalpart1 + fac*log(betasig) + gammaln(0.5*nstar) - gammaln(alphasig);
    
    if inc_int == 1
        bigX = obj.X1;
        if size(obj.X1,2) > 1
            bigX = cat(1, zeros(obj.n1, size(obj.X0,2)),obj.X0);
            bigX = cat(2, bigX, obj.X1);
        end
        Hvalpart1 = - fac*log(2*pi) + fac*log(betasig) + gammaln(0.5*nstar) ...
            - gammaln(alphasig) - 0.5*log(det(bigX'*bigX));
    end
end
end

function [Hvalpart2] = getHpart2(NpointinTimeIn, integrated, inc_int, betasig)
obj = computefromInt(NpointinTimeIn);
Hvalpart2 = obj.Sigma - obj.Sigma*obj.X1*inv(obj.X1'*obj.Sigma*obj.X1)*obj.X1'*obj.Sigma;
if integrated == 1
    Hvalpart2 = Hvalpart2 * betasig * 0.5;
    if inc_int == 1
        bigX = obj.X1;
        if size(obj.X1,2) > 1
            bigX = cat(1, zeros(obj.n1, size(obj.X0,2)),obj.X0);
            bigX = cat(2, bigX, obj.X1);
        end
        Hvalpart2 = betasig * 0.5 * (eye(size(bigX,1)) - bigX*inv(bigX'*bigX)*bigX');
    end
end
end

function [Hval] = getHpartInfo(NpointinTimeIn)
global inc_ints inc_int
obj = computefromInt(NpointinTimeIn);
Hval.Sigma = obj.Sigma;
Hval.logdetX0 = log(det(obj.X0'*obj.X0));
Hval.X1Sigma = obj.X1'*obj.Sigma;
Hval.invX1SigmaX1 = inv(Hval.X1Sigma*obj.X1);
Hval.logdetX1SigmaX1 = -log(det(Hval.invX1SigmaX1));
if inc_int == 1
    Hval.bigX = obj.X1;
    if size(obj.X1,2) > 1
        Hval.bigX = cat(1, zeros(obj.n1, size(obj.X0,2)),obj.X0);
        Hval.bigX = cat(2, Hval.bigX, obj.X1);
    end
    Hval.XX = inv(Hval.bigX' * Hval.bigX);
end
if inc_ints==1
    bigX = [obj.X0, obj.X1];
    Hval.bigXX = inv(bigX' * bigX);
    Hval.XX = inv(obj.X0'*obj.X0);
end
end

function [x] = UpdateAlpha0(x)
global sitesize aAlpha0 bAlpha0

alpha = x{1}.alpha;
Nstar = size(x,1);
eta = betarnd(alpha+1,sitesize);
pieta = (aAlpha0+Nstar-1)/(aAlpha0+Nstar-1+sitesize*(1/bAlpha0-log(eta)));
u = rand(1);
if u > pieta
    x{1}.alpha = gamrnd(aAlpha0+Nstar, 1/(1/bAlpha0 - log(eta)));
else
    x{1}.alpha = gamrnd(aAlpha0+Nstar-1, 1/(1/bAlpha0 - log(eta)));
end
end

function [x] = UpdateClusterInfoNofixedsigma(x)
global MinNtimesw timesize Nint inc_int inc_ints Hinfo

Nk = Nint;
runtime = cputime;
Ncluster = size(x,1);

for nocluster = 1:Ncluster
    B = cell(Nk,1);
    vbetaltemp = cell(Nk,1);
    mubeta1temp = cell(Nk,1);
    SumH = cell(Nk,1);
    SumHPluslogmultinomial = cell(Nk,1);

    sizedata = x{nocluster,1}.Siteincluster;
    sizeN = size(sizedata,2);
    tmpV = [];
    %MatV=[];
    for i = 1:sizeN
        tmpV = cat(1,tmpV,(x{nocluster,1}.data{i,1}.TimeIntervals)');
    end
    
    MatV= cell2mat(tmpV);% all data value in the same cluster (stored by matrix form)
    
    %B=[];SumH=[];Sumdata=[];
    for Ni=1:Nk
        B{Ni,1} = Combination(Ni,timesize-Ni*MinNtimesw);%all posssible combination
        m = size(B{Ni,1},1); k = size(B{Ni,1},2);
        for s = 1:m,
            [t1,t2,t3] = Hhat(x{nocluster,1},B{Ni,1},s,sizeN);%for sdd=1:sizeN,
            SumH{Ni,1}(s) = t1; vbetaltemp{Ni,1}{s} = t2; mubeta1temp{Ni,1}(s,:) = t3;
        end
        newB = B{Ni,1}-MinNtimesw;
        SumHPluslogmultinomial{Ni,1} = SumH{Ni,1}'+logmultinomial(newB,ones(1,k)./Ni);
    end
    %generate number of time interval
    
    UNlogProbK = cellfun(@LogMAxfunction,SumHPluslogmultinomial)'+logtruncpoissonpdf(((1:Nk)-1),x{nocluster,1}.lambda);
    ProbK = exp(UNlogProbK-max(UNlogProbK))./sum(exp(UNlogProbK-max(UNlogProbK)));
    IndexK = randsample(1:Nk,1,true,ProbK);
    
    %generate how many data points in each time interval
    
    UNlogProbNK = SumHPluslogmultinomial{IndexK,1};
    ProbNK = exp(UNlogProbNK-max(UNlogProbNK))./sum(exp(UNlogProbNK-max(UNlogProbNK)));
    
    Index = randsample(1:numel(ProbNK),1,true,ProbNK);
    
    nV = B{IndexK,1}(Index,:);
    for j = 1:sizeN,
        x{nocluster,1}.data{j,1}.TimeIntervals = mat2cell(MatV(j,:),1,nV)';
    end
    x{nocluster,1}.NpointinTimeIn = nV;
    x{nocluster,1}.NTimeIntervals = numel(nV);
    %generate theta in each time interval
    %     randn('state',bloop); bloop = bloop+10;
    if inc_int ~= 1
        x{nocluster,1}.parameter = mvnrnd(mubeta1temp{IndexK,1}(Index,:),vbetaltemp{IndexK,1}{Index});
        
        if inc_ints == 1
            y = x{nocluster,1}; obj = computefromInt(nV);
            for sdd = 1:sizeN
                data = cat(2,y.data{sdd}.TimeIntervals{:,1}) - y.u(sdd) - (obj.X1*x{nocluster,1}.parameter')';
                vbeta0temp = y.sigma2(sdd).*Hinfo{IndexK}{Index}.XX;
                mubeta0temp = vbeta0temp*obj.X0'*(data')./y.sigma2(sdd);
                x{nocluster,1}.beta0{sdd} = mvnrnd(mubeta0temp,vbeta0temp);
            end
        end
    else
        betas = mvnrnd(mubeta1temp{IndexK,1}(Index,:),vbetaltemp{IndexK,1}{Index});
        x{nocluster,1}.beta0 = [0 betas(1:(IndexK-1))];
        x{nocluster,1}.parameter = betas(IndexK:end);
   end
end
runtime = cputime-runtime;
end

function [x] = UpdateGamma(x)
global W u0 loglike0 gap gammas
time = cputime;
A = stackU(x) - u0';
tmp = A'*W*A./(2*x{1}.tau2); 
loglike = loglike0 + gammas*tmp; 
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
u = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(u > cump);
gammanew = gammas(1);
if(i0>1)
    gammanew = gammas(i0-1) + gap/P(i0)*(u-cump(i0));
end

Nstar = size(x,1);
cellgammaNew = mat2cell(gammanew*ones(Nstar,1), ones(Nstar,1));
[x]=cellfun(@replacegamma,x,cellgammaNew,'UniformOutput',false);

runtime=cputime-time;
end

function [x] = replacegamma(x,gammaNew0)
x.gamma = gammaNew0;
end

function [x] = UpdatelambdaNofixedsigma(x)
global alphalambda  betalambda kstar sitesize fixsmalllambda
time = cputime;
k = cellfun(@CountTimeInterval,x);
k = k-1;
Nstar = size(k,1);
Nsitesincluster = cellfun(@Countsitenumber,x);
sumk = k'*Nsitesincluster;

gap = 1e-2; 
% previous: 1e-4:1e-2:1e-1
lambda0 = x{1}.lambda;
if fixsmalllambda == 0
    lambdas = max(gap,lambda0-200*gap):gap:(lambda0+200*gap);
else
    lambdas = 1e-4:1e-2:1e-1;
end

lvec = 0:kstar;
loglike = (alphalambda-1+sumk)*log(lambdas) -lambdas/betalambda - sitesize*log(truncPoi(lambdas, lvec));
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
% plot(lambdas, P);
u = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(u > cump);
lambdaNew = lambdas(1);
if(i0>1)
    lambdaNew = lambdas(i0-1) + gap/P(i0)*(u-cump(i0));
end

celllambdaNew = mat2cell(lambdaNew*ones(Nstar,1), ones(Nstar,1));
[x]=cellfun(@replacelambda,x,celllambdaNew,'UniformOutput',false);

runtime=cputime-time;
end 

function [tsum] = truncPoi(lambdas, lvec)
tsum = zeros(1,length(lambdas));
for i = 1:length(lambdas)
  tsum(i) = sum(lambdas(i).^lvec./factorial(lvec));
end
end

function [x] = replacelambda(x,lambdaNew0)
x.lambda=lambdaNew0;
end

function [logq0val,Bk,SumHaddLogMulnomialLogpoisson] = Updateq0(SumHaddLogMulnomial,lambda)
global sitesize Nint

SumHaddLogMulnomialLogpoisson = cell(Nint,sitesize);
for j=1:sitesize,
    for i=1:Nint % the number of time intervals
        SumHaddLogMulnomialLogpoisson{i,j} = SumHaddLogMulnomial{i,j} + logtruncpoissonpdf(i-1,lambda);
    end
end
Bk = cellfun(@LogMAxfunction,SumHaddLogMulnomialLogpoisson);
%q0=sum(exp(Bk),1);

logq0val = zeros(1,sitesize);
for j=1:sitesize
    logq0val(j) = LogMAxfunction(Bk(:,j));
end
%NormalizeBk=exp(Bk-ones(N,1)*max(Bk))./(ones(N,1)*sum(exp(Bk-ones(N,1)*max(Bk)),1));
end

function [x] = Updatesigma2Nofixedsigma(x)
global sitesize timesize alphasig betasig
global priorchoice inc_int inc_ints

time = cputime;
Nstar = size(x,1); 
Sigmas = cell(1,Nstar); X1s = cell(1,Nstar); Kss = zeros(1,Nstar);
X0s = cell(1,Nstar); N1s = cell(1,Nstar);

clLabel = []; marker1 = []; marker2 = [];
for i = 1:Nstar
    tmp = computefromInt(x{i}.NpointinTimeIn);
    Sigmas{i} = tmp.Sigma;
    X1s{i} = tmp.X1;
    if inc_int==1
        X0s{i} = tmp.X0; N1s{i} = tmp.n1;
    elseif inc_ints==1
        X0s{i} = tmp.X0;
    end
    Kss(i) = length(x{i}.NpointinTimeIn) - 1;
    clLabel = [clLabel, x{i}.Siteincluster];
    marker1 = [marker1, i*ones(1,length(x{i}.Siteincluster))];
    marker2 = [marker2, 1:length(x{i}.Siteincluster)];
end
[clLabel, I] = sort(clLabel); marker1 = marker1(I); marker2 = marker2(I);

for j=1:sitesize
    c = marker1(j); yy = marker2(j);
    if priorchoice == 1
        vec = cat(2,x{c}.data{yy}.TimeIntervals{:,1})' - x{c}.u(yy) - X1s{c}*x{c}.parameter';
        a0 = timesize - Kss(marker1(j));
        azeros = 0.5*a0 + alphasig;
        bzeros = (0.5*vec'*Sigmas{c}*vec + betasig^-1).^-1;
        
        if inc_int==1
            vec((N1s{c}+1):end) = vec((N1s{c}+1):end) - X0s{c}*x{c}.beta0(2:end)';
            azeros = 0.5*timesize + alphasig;
            bzeros = (0.5*(vec'*vec) + betasig^-1).^-1;
        elseif inc_ints==1
            vec = vec - X0s{c}*x{c}.beta0{yy}';
            azeros = 0.5*timesize + alphasig;
            bzeros = (0.5*(vec'*vec) + betasig^-1).^-1;
        end
        
        x{c,1}.sigma2(yy) = 1./gamrnd(azeros,bzeros);
    end
end
end

function [x] = UpdatesiteNofixedsigma(x)
global sitesize MinNtimesw timesize  SumHpart1addLogMulnomial 
global Hpart2 Nint Hinfo integrated alphasig betasig inc_int inc_ints Lungdata
% global marker1 marker2 obj

time = cputime;
alpha = x{1,1}.alpha;
lambda = x{1,1}.lambda;

Nstar = size(x,1); 
clLabel = []; marker1 = []; marker2 = [];
for i = 1:Nstar
    clLabel = [clLabel, x{i}.Siteincluster];
    marker1 = [marker1, i*ones(1,length(x{i}.Siteincluster))];
    marker2 = [marker2, 1:length(x{i}.Siteincluster)];
end
[clLabel, I] = sort(clLabel); marker1 = marker1(I); marker2 = marker2(I); 

% combine SumHpart1addLogMulnomial and Hpart2 to SumHaddLogMulnomial
SumHaddLogMulnomial = cell(Nint, sitesize);
for j = 1:sitesize
    tmpsigma2 = x{marker1(j)}.sigma2(marker2(j));
    tmpY = cat(2,x{marker1(j)}.data{marker2(j)}.TimeIntervals{:,1}) - x{marker1(j)}.u(marker2(j));
    for i = 1:Nint  %note i = Ks + 1
        tmps = [];
        for l = 1:length(Hpart2{i}) % number of combinations for this number of intervals
            if integrated == 0
                SumHaddLogMulnomial{i,j}(l,1) = SumHpart1addLogMulnomial{i}(l,1) ...
                    - 0.5*(timesize-2*i+1)*log(tmpsigma2) ...
                    - 0.5/tmpsigma2*tmpY*Hpart2{i}{l}*tmpY';
            elseif integrated == 1
                tmps(l) = tmpY*Hpart2{i}{l}*tmpY';
                SumHaddLogMulnomial{i,j}(l,1) = SumHpart1addLogMulnomial{i}(l,1) ...
                    - 0.5*(2*alphasig+timesize-2*i+1)*log(1+tmps(l));
            end
        end
    end
end

[logq0,Bk,SumHaddLogMulnomialLogpoisson] = Updateq0(SumHaddLogMulnomial,lambda);

for j = 1:sitesize %sitesize
    clLabel = []; marker1 = []; marker2 = []; % need to update label information
    Nstar = size(x,1);% define N*: the number of cluster
    N = zeros(1,Nstar);
    logdets = cell(1,Nstar); Xbetas = cell(1,Nstar); Sigmas = cell(1,Nstar);
    Logq0 = log(alpha)+logq0(j);
    for i = 1:Nstar
        N(i) = size(x{i}.Siteincluster,2);% N(i): the number of sites in cluster i.
        tmp = computefromInt(x{i}.NpointinTimeIn);
        logdets{i} = -0.5*log(det(tmp.X0'*tmp.X0));
        Xbetas{i} = tmp.X1*x{i}.parameter';
        if inc_int == 1
            Xbetas{i}((tmp.n1+1):end) = Xbetas{i}((tmp.n1+1):end) + tmp.X0*x{i}.beta0(2:end)'; 
        end
        Sigmas{i} = tmp.Sigma;
        clLabel = [clLabel, x{i}.Siteincluster];
        marker1 = [marker1, i*ones(1,length(x{i}.Siteincluster))];
        marker2 = [marker2, 1:length(x{i}.Siteincluster)];
    end
    [clLabel, I] = sort(clLabel); marker1 = marker1(I); marker2 = marker2(I);
    c = marker1(j); yy = marker2(j);
    tmpsigma2 = x{c}.sigma2(yy);
    data = cat(2,x{c}.data{yy}.TimeIntervals{:,1}) - x{c}.u(yy);
    
    qp = zeros(1,Nstar);
    for i = 1:Nstar
        if integrated == 0
            qp(i) = -0.5*(timesize - x{i}.NTimeIntervals) * log(2*pi*tmpsigma2) + logdets{i}...
                - 0.5/tmpsigma2*(data-Xbetas{i}')*Sigmas{i}*(data'-Xbetas{i});
        elseif integrated == 1
            fac = 0.5 * (timesize - x{i}.NTimeIntervals);
            if inc_int == 1
                fac = 0.5 * timesize; logdets{i} = 0; Sigmas{i} = 1;
            end
            qp(i) = - fac*log(2*pi) + logdets{i} + gammaln(fac + alphasig) - gammaln(alphasig) - alphasig*log(betasig)...
                - (fac + alphasig)*log(0.5*(data-Xbetas{i}')*Sigmas{i}*(data'-Xbetas{i}) + 1/betasig);
        end
    end
    N(c) = N(c)-1;
    MaxQ = max([Logq0 qp+log(N)]);
    SumQ = sum(exp([Logq0 qp+log(N)]-MaxQ));
    qprob = exp([Logq0 qp+log(N)]-MaxQ)/SumQ;
    Newlevel = randsample(0:Nstar,1,true,qprob);
    if (Newlevel==c),% j site does not move
        N(c)=N(c)+1;
    else
        % move j from the cluster c to the others 
        if (Newlevel==0),%Create new cluster and j will be moved to the new cluster
            Nstar = Nstar+1;
            NTimeIntProb = exp(Bk(:,j)-max(Bk(:,j)))./sum(exp(Bk(:,j)-max(Bk(:,j))));
            NewNTimeInt = randsample(1:size(Bk(:,j),1),1,true,NTimeIntProb);
            %generate how many time points in each interval
            NewB = Combination(NewNTimeInt,timesize-NewNTimeInt*MinNtimesw);
            UNlogProbInd = SumHaddLogMulnomialLogpoisson{NewNTimeInt,j};
            logProbInd = exp(UNlogProbInd-max(UNlogProbInd))./sum(exp(UNlogProbInd-max(UNlogProbInd)));
            NewTimeInd = randsample(1:size(NewB,1),1,true,logProbInd);
            
            x{Nstar,1}.NpointinTimeIn = NewB(NewTimeInd,:);
            
            if inc_int ~= 1
                if inc_ints ~= 1
                    newvlsite = Hinfo{NewNTimeInt}{NewTimeInd}.invX1SigmaX1*tmpsigma2;
                    newbetasite = newvlsite*Hinfo{NewNTimeInt}{NewTimeInd}.X1Sigma*data'/tmpsigma2;
                    x{Nstar,1}.parameter = mvnrnd(newbetasite, newvlsite);
                else
                    tmp = computefromInt(x{Nstar,1}.NpointinTimeIn); tmp = [tmp.X0, tmp.X1];
                    newvlsite = Hinfo{NewNTimeInt}{NewTimeInd}.bigXX*tmpsigma2;
                    newbetasite = newvlsite*tmp'*data'/tmpsigma2;
                    tmpeBetas = mvnrnd(newbetasite, newvlsite);
                    x{Nstar,1}.beta0{1} = tmpeBetas(1:NewNTimeInt);
                    x{Nstar,1}.parameter = tmpeBetas((NewNTimeInt+1):end);
                end
            else
                newvlsite = Hinfo{NewNTimeInt}{NewTimeInd}.XX*tmpsigma2;
                newbetasite = newvlsite*Hinfo{NewNTimeInt}{NewTimeInd}.bigX'*data'/tmpsigma2;
                tmpeBetas = mvnrnd(newbetasite, newvlsite);
                x{Nstar,1}.beta0 = [0 tmpeBetas(1:(NewNTimeInt-1))];
                x{Nstar,1}.parameter = tmpeBetas(NewNTimeInt:end);
            end
            x{Nstar,1}.Siteincluster = j;
            x{Nstar,1}.sigma2 = x{c,1}.sigma2(yy); %directly copy
            x{Nstar,1}.NTimeIntervals = NewNTimeInt; 
            x{Nstar,1}.data{1,1}.TimeIntervals = mat2cell(data,1,x{Nstar,1}.NpointinTimeIn)';
            x{Nstar,1}.u = x{c,1}.u(yy); %directly copy
            x{Nstar,1}.lambda = x{c,1}.lambda;
            x{Nstar,1}.tau2 = x{c,1}.tau2;
            x{Nstar,1}.gamma = x{c,1}.gamma;
        else %j will be moved to the existed cluster.
            x{Newlevel,1}.Siteincluster=[x{Newlevel,1}.Siteincluster j];
            N(Newlevel)=N(Newlevel)+1;
            x{Newlevel,1}.data{N(Newlevel),1}.TimeIntervals= mat2cell(data,1,x{Newlevel,1}.NpointinTimeIn)';
            x{Newlevel,1}.sigma2(N(Newlevel)) = x{c,1}.sigma2(yy); %directly copy
            x{Newlevel,1}.u(N(Newlevel)) = x{c,1}.u(yy); %directly copy
            if inc_ints == 1
                obj = computefromInt(x{Newlevel,1}.NpointinTimeIn);
                vbeta0temp = tmpsigma2.*((obj.X0'*obj.X0)\eye(size(obj.X0,2)));
                mubeta0temp = vbeta0temp*obj.X0'*(data'-obj.X1*x{Newlevel,1}.parameter')./tmpsigma2;
                x{Newlevel,1}.beta0{N(Newlevel)} = mvnrnd(mubeta0temp,vbeta0temp);
            end
        end

        % Cancel the site j's information in Cluster c after it is moved to the other cluster.
        x{c,1}.Siteincluster(yy)=[];
        x{c,1}.sigma2(yy)=[];
        x{c,1}.u(yy)=[];
        x{c,1}.beta0=[x{c,1}.beta0(1:(yy-1)), x{c,1}.beta0((yy+1):end)];
        % after moving, if the cluster c become empty, then cancel the cluster c.
        if (size(x{c,1}.Siteincluster,2)==0),
            if (c==1),
                x=x(2:end,:);
            elseif (c==Nstar),
                x=x(1:end-1,:);
            else
                x=[x(1:c-1,:) ;x(c+1:end,:)];
            end
            Nstar = Nstar-1;
            N(c)=[];
        else
            if (yy==1),
                x{c,1}.data=x{c,1}.data(2:end,:);
            elseif (yy==size(x{c,1}.data,1)),
                x{c,1}.data=x{c,1}.data(1:end-1,:);
            else
                x{c,1}.data=[x{c,1}.data(1:yy-1,:) ;x{c,1}.data(yy+1:end,:)];
            end
        end
    end
end

x{1,1}.alpha = alpha;
end

function [x] = UpdateTau2(x)
global alphatau2 betatau2 W M sitesize u0
% global w0 Sigma2Range

time = cputime;
t = stackU(x);
V = (M - x{1}.gamma*W);
A = 0.5*(t-u0')'*V*(t-u0');

sigmaprior = 2;
switch sigmaprior
    case 1 %default prior   % NOT YET MODIFIED
        
        sigmavec = Sigma2Range;
        N = length(sigmavec);
        loglike = zeros(1,N);
        for i = 1:N
            sig2 = sigmavec(i);
            logprior=-2*log(1+w0*sig2);
            loglike(i) = -A/sig2 - (n/2)*log(sig2) +logprior;
        end
        
        MaxLogLike = max(loglike);
        P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
        %plot(sigmavec,P,'k-')
        %plot(sigmavec,loglike,'k-')
        u = rand(1);
        bloop = bLoc;
        cump = cumsum([0 P(1:(end-1))]);
        i0 = sum(u > cump);
        X.P{3} = sigmavec(1);
        if(i0>1)
            X.P{3} = sigmavec(i0-1) + .001/P(i0)*(u-cump(i0));
        end
        %X.P{3} = sigmavec(i0);
        
    case 2 % gamma prior
        a = sitesize/2+alphatau2;
        b = 1/(A + 1/betatau2);
        y = gamrnd(a,b);
        tau2new  = 1/y;
        Nstar = size(x,1);
        celltau2New = mat2cell(tau2new*ones(Nstar,1), ones(Nstar,1));
        [x] = cellfun(@replacetau2,x,celltau2New,'UniformOutput',false);
end
end

function [x] = replacetau2(x,tau2New0)
x.tau2 = tau2New0;
end

function [x] = UpdateU(x)
global M W sitesize timesize u0
global priorchoice delta2 inc_int inc_ints

reducedX0 = 1;

time = cputime;
invtauD = 1/x{1}.tau2 * (M - x{1}.gamma*W);
A = zeros(sitesize);
B = zeros(sitesize, 1);

Nstar = size(x,1);
obj = cell(1,Nstar); clLabel = []; marker1 = []; marker2 = [];
for i = 1:Nstar
    obj{i} = computefromInt(x{i}.NpointinTimeIn);
    clLabel = [clLabel, x{i}.Siteincluster];
    marker1 = [marker1, i*ones(1,length(x{i}.Siteincluster))];
    marker2 = [marker2, 1:length(x{i}.Siteincluster)];
end
[clLabel, I] = sort(clLabel); marker1 = marker1(I); marker2 = marker2(I);


for i = 1:sitesize
    tmpsigma2 = x{marker1(i)}.sigma2(marker2(i));
    tmpY = cat(2,x{marker1(i)}.data{marker2(i)}.TimeIntervals{:,1});
    if(~reducedX0)
        if priorchoice == 1
            tmpSigma = obj{marker1(i)}.Sigma;
        elseif priorchoice == 2
            tmpX0 = obj{marker1(i)}.X0;
            tmpSigma = eye(timesize) - tmpX0*inv(tmpX0'*tmpX0 + 1/delta2*tmpsigma2*eye(size(tmpX0,2)))*tmpX0';
        end
        A(i,i) = ones(1,timesize)*tmpSigma*ones(timesize,1)/tmpsigma2;
        B(i,1) = ones(1,timesize)*tmpSigma*(tmpY'-obj{marker1(i)}.X1*...
            x{marker1(i)}.parameter')/tmpsigma2;
    elseif(reducedX0 == 1)
        tmp = tmpY'-obj{marker1(i)}.X1*x{marker1(i)}.parameter';
        n1 = obj{marker1(i)}.n1;
        if inc_int == 1
            tmp((n1+1):end) = tmp((n1+1):end) - obj{marker1(i)}.X0*x{marker1(i)}.beta0(2:end)';
            n1 = timesize;
        elseif inc_ints == 1
            tmp = tmp - obj{marker1(i)}.X0*x{marker1(i)}.beta0{marker2(i)}';
            n1 = timesize;
        end
        A(i,i) = n1/tmpsigma2;
        B(i,1) = A(i,i) * sum(tmp(1:n1))/n1;
    end
end

SigmaU = (A + invtauD)\eye(size(A));
muU = SigmaU * (B + invtauD*u0');
Unew = mvnrnd(muU, SigmaU);
for i = 1:sitesize
    x{marker1(i)}.u(marker2(i)) = Unew(i);
end
runtime=cputime-time;
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

function [k] = CountTimeInterval(x)
k=x.NTimeIntervals;
end

function [SumHval,vbetalnew,mubeta1new] = Hhat(y,Ballcombination,combinationN,NSite)
global Hinfo timesize inc_int inc_ints

NtimeInt = size(Ballcombination, 2);

SumHval = - 0.5*(timesize-NtimeInt*NSite - NtimeInt)*log(2*pi) ...
    - 0.5*(timesize-NtimeInt+1)*sum(log(y.sigma2))...
    - 0.5*NSite*Hinfo{NtimeInt}{combinationN}.logdetX0...
    - 0.5*NtimeInt*log(sum(y.sigma2.^(-1))) - 0.5*Hinfo{NtimeInt}{combinationN}.logdetX1SigmaX1;

% individual information
sum1 = 0.0; sum2 = 0.0; sum15 = 0.0;
Sigma = Hinfo{NtimeInt}{combinationN}.Sigma;
X1Sigma = Hinfo{NtimeInt}{combinationN}.X1Sigma;
invX1SigmaX1 = Hinfo{NtimeInt}{combinationN}.invX1SigmaX1;

for sdd = 1:NSite,
    data = cat(2,y.data{sdd}.TimeIntervals{:,1}) - y.u(sdd);
    sum1 = sum1 + data*Sigma*data'/y.sigma2(sdd);
    sum15 = sum15 + data'/y.sigma2(sdd);
    sum2 = sum2 + X1Sigma * data'/y.sigma2(sdd);
end
SumHval = SumHval - 0.5*sum1 + 0.5 *sum2'*invX1SigmaX1*sum2/sum(y.sigma2.^(-1));

if inc_int ~= 1
    vbetalnew = invX1SigmaX1/sum(y.sigma2.^(-1));
    mubeta1new = vbetalnew*sum2;
else
    vbetalnew = Hinfo{NtimeInt}{combinationN}.XX/sum(y.sigma2.^(-1));
    mubeta1new = vbetalnew*Hinfo{NtimeInt}{combinationN}.bigX'*sum15;
end
end

function [y] = LogMAxfunction(x)
% compute log(sum(exp(x)))
y=log(sum(exp(x-max(x))))+max(x);
end

function [y] = logmultinomial(x,y)
y=gammaln(sum(x,2)+1)-sum(gammaln(x+1),2)+x*log(y)';
end

function [Val] = logtruncpoissonpdf(x,lambda)
global kstar
%Function to calculate the log of iid truncated poisson pdf
vec = 0:kstar;
Val = x*log(lambda) - gammaln(x+1) - log(sum(lambda.^vec./factorial(vec)));
end

function [U] = stackU(x)
global sitesize
clLabel = []; marker1 = []; marker2 = [];
Nstar = size(x,1);
for i = 1:Nstar
    clLabel = [clLabel, x{i}.Siteincluster];
    marker1 = [marker1, i*ones(1,length(x{i}.Siteincluster))];
    marker2 = [marker2, 1:length(x{i}.Siteincluster)];
end
[clLabel, I] = sort(clLabel); marker1 = marker1(I); marker2 = marker2(I); 
U = zeros(sitesize,1); 
for i = 1:sitesize
    U(i,1) = x{marker1(i)}.u(marker2(i));
end
end

function [out] = RessCreatedata(datavalue,x)
[betas, bint, Ressval] =regress(datavalue,x);
out.betazero = betas(1);
out.betaone = betas(2); 
out.Ressval = Ressval;
end

function [ksites] = Countsitenumber(x)
ksites=size(x.Siteincluster,2);
end
