function [] = outH2Rdata()
global sitesize
dirs = './';
rep = 1; c = 6;
datachoice = 1;
load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat'))   % get true_bvec
data = data(2:end,:);
t = cat(1,1:size(data,2),data);
timesize = size(data,1);
load(strcat(dirs,'outH_',num2str(rep),'.mat'));
Dist = outH.Dist; 
sitesize = size(Dist,1);
pdistvec=[];
for i=1:sitesize-1,
    pdistvec = [pdistvec,Dist(i,(i+1):end)];
end
Dmin = min(pdistvec); Dmax = max(pdistvec);
pdistvec = (pdistvec-Dmin)/(Dmax-Dmin);
clusterRe = linkage(pdistvec, 'ward');
% Name=char('Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Washington DC','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky',...
%     'Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire	',...
%     'New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina',...
%     'South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming');
% [H,T] = dendrogram(clusterRe,sitesize, 'LABELS', Name,'orient','left','colorthreshold',1.3);
% set(H,'LineWidth',2);axis tight;
% title('')
% set(gca,'FontSize',13)
% orient landscape
% print('-painters', '-dpsc2', '-r600', strcat('DendrogramP_',num2str(rep),'.eps'))

ClusterIndex = cluster(clusterRe,'maxclust',c);
Clustergroup = cell(c,1);
for i = 1:sitesize,
    for j = 1:c
        if  ClusterIndex(i)==j,
            Clustergroup{j,1} = [Clustergroup{j,1} i];
        end
    end
end
clusterMeasure = computefromTable(Clustergroup, Clustergroup0);
disp(clusterMeasure.tab)

Rdatas = cell(1,c);
for r = 1:c,
    allmonitorsite = Clustergroup{r,1}; nr = size(allmonitorsite,2);
    Rdata.sites = allmonitorsite;
    Rdata.ObsY = zeros(timesize,nr);
    Rdata.PredY = zeros(timesize,nr);
    Rdata.miny = zeros(1,nr);
    Rdata.maxy = zeros(1,nr);
    Rdata.sigma2 = zeros(1,nr);
    Rdata.yvalueupperCI = zeros(timesize,nr);
    Rdata.yvaluelowerCI = zeros(timesize,nr);
    Rdata.cpt = cell(1,nr); 
    for monidex = 1:size(allmonitorsite,2)
        yhat = [];
        Rdata.ObsY(:,monidex) = t(2:end,allmonitorsite(monidex));
        B = outH.Hcps{allmonitorsite(monidex)};
        TimeInterval = [0,cumsum(B)];
        Rdata.cpt{monidex} = zeros(1,size(TimeInterval,2) - 1);
        for i = 1:size(TimeInterval,2) - 1,
            Rdata.cpt{monidex}(i) = TimeInterval(i)+1; vec = (TimeInterval(i)+1):TimeInterval(i+1);
            yhat = [yhat, outH.Hbeta0s{allmonitorsite(monidex)}(i) + vec*outH.Hbeta1s{allmonitorsite(monidex)}(i)];
        end
        Rdata.PredY(:,monidex) = yhat';
        Rdata.miny(monidex) = min([Rdata.ObsY(:,monidex)', Rdata.PredY(:,monidex)']);
        Rdata.maxy(monidex) = max([Rdata.ObsY(:,monidex)', Rdata.PredY(:,monidex)']);
    end
    Rdatas{r} = Rdata;
end
save(strcat('RdataP',num2str(rep),'.mat'),'Rdatas')
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
