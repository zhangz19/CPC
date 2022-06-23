function [] = getTables()

% for simulation table
cases = 1:12;
model = 3; % 49 states
dataPath = './'; %simu3/';
for rep = cases
    load(strcat(dataPath,'PsOutputModel',num2str(model),'_',num2str(rep),'.mat'))
    fprintf('$%.d$',rep); 
    paras = [mean(cptMeasure.accuracy), clusterMeasure.SS, clusterMeasure.SS2,...
        mean(outH.HRa), outH.SS2(3,4), outH.SS2(4,4), outH.SS2(5,4), outH2.SS2(3,4), outH2.SS2(4,4), outH2.SS2(5,4)];
    fprintf(repmat('&$%.3f$', [1,length(paras)]), paras)
    fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
end
