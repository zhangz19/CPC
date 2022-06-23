function [] = getMeasureSimple()
% explore the clustering results from penalized methods with simple feature
% matrix formed by replicating beta over time, to address one reviewer's concern
% add outH2 to PsOutputModel
cases = 12;
simu = 1; datachoice = 1;
model = 3; % 49 states
dataPath = './';
timesize = 38; sitesize = 49;

for rep = 1:cases
    disp(rep)
    D = nan(sitesize, timesize);
    
    load(strcat('simulateddata',num2str(datachoice),'_',num2str(rep),'.mat')) %Clustergroup0
    load(strcat(dataPath,'PsOutputModel',num2str(model),'_',num2str(rep),'.mat'))
    for s = 1:sitesize
        tmp0 = cumsum([1 outH.Hcps{s}]);
        for l = 1:(length(tmp0)-1)
            D(s,tmp0(l):(tmp0(l+1)-1)) = outH.Hbeta1s{s}(l);
        end
    end
    
    D1 = zeros(sitesize);
    for s1 = 2:sitesize % find clustering measures
        for s2 = 1:(s1-1)
            D1(s1,s2) = sqrt(sum( (D(s1,:)-D(s2,:)).^2 ));
        end
    end
    for s1 = 2:sitesize
        for s2 = 1:(s1-1)
            D1(s2,s1) = D1(s1,s2);
        end
    end
    pdistvec=[];
    for i=1:sitesize-1,
        pdistvec = [pdistvec D1(i,i+1:end)];
    end
    clusterRe = linkage(pdistvec, 'ward');
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
        end
        out.SS2 = [cvec', Sensitivity, Specificity, 2-Sensitivity-Specificity];
    end    
    outH2 = out; 
    save(strcat(dataPath,'PsOutputModel',num2str(model),'_',num2str(rep),'.mat'), ...
        'PsOutput','fitMeasure','HPD','R','clusterMeasure','cptMeasure','Clustergroup','outH', 'outH2')
end