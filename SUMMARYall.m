function [results] = SUMMARYall()

cases = 6;
REP = 1; M = 3; sitesize = 49;  % 49 states
dataPath = './'; simu=1;

% sitewise results
results.fit.DIC = cell(1,sitesize);
results.fit.BIC = cell(1,sitesize);
results.fit.SSE = cell(1,sitesize);
results.fit.MSE = cell(1,sitesize);
results.cpt.HPD = cell(1,sitesize);
results.cpt.KL = cell(1,sitesize); % KL divergence between posterior distribution of changepoint and the true mass
results.cpt.entropy = cell(1,sitesize); % entropy of posterior
results.cpt.relentropy = cell(1,sitesize); % relative entropy of posterior
results.cpt.accuracy = cell(1,sitesize); % mean accuracy rate for detection of true changepoints
results.cpt.betadiff = cell(1,sitesize); % difference of beta
results.cpt.betabwHPD = cell(1,sitesize); % HPD band width of difference
results.cpt.betabwCI = cell(1,sitesize); % CI band width of difference

for model = 3:M
    
    for rep = 1:cases
        load(strcat(dataPath,'PsOutputModel',num2str(model),'_',num2str(rep),'.mat')) % cases instead of rep
        
        if simu == 1
            results.cl.chisq(rep, model) = clusterMeasure.chisq;
            results.cl.phi(rep, model) = clusterMeasure.phi;
            results.cl.yuleQ(rep, model) = clusterMeasure.yuleQ;
            results.cl.specificity(rep, model) = clusterMeasure.specificity;
            results.cl.sensitivity(rep, model) = clusterMeasure.sensitivity;
            results.cl.Cohenkappa(rep, model) = clusterMeasure.Cohenkappa;
            results.cl.SS2(rep, model) = clusterMeasure.SS2;
        end
        results.cl.entropy(rep, model) = clusterMeasure.entropy;
        results.cl.relentropy(rep, model) = clusterMeasure.relentropy;
        results.cl.SS(rep, model) = clusterMeasure.SS;
        
        for s = 1:sitesize
            results.fit.DIC{s}(rep, model) = fitMeasure.DIC(s);
            results.fit.BIC{s}(rep, model) = fitMeasure.BIC(s);
            results.fit.SSE{s}(rep, model) = fitMeasure.SSE(s);
            results.fit.MSE{s}(rep, model) = fitMeasure.MSE(s);
            results.cpt.HPD{s}{rep, model} = cptMeasure.HPD{s};
            results.cpt.entropy{s}(rep, model) = cptMeasure.entropy(s);
            results.cpt.relentropy{s}(rep, model) = cptMeasure.relentropy(s);
            if simu == 1
                results.cpt.KL{s}(rep, model) = cptMeasure.KL(s);
                results.cpt.accuracy{s}(rep, model) = cptMeasure.accuracy(s);
            end
            results.cpt.betadiff{s}(rep, model) = cptMeasure.betadiff(s);
            results.cpt.betabwHPD{s}(rep, model) = cptMeasure.betabwHPD(s);
            results.cpt.betabwCI{s}(rep, model) = cptMeasure.betabwCI(s);
        end
    end
    
%     out.mDIC(model) = mean(results.DIC(:, model));
%     out.mcl.chisq(model) = mean(results.cl.chisq(:, model));
%     out.mcl.phi(model) = mean(results.cl.phi(:, model));
%     out.mcl.yuleQ(model) = mean(results.cl.yuleQ(:, model));
%     out.mcl.specificity(model) = mean(results.cl.specificity(:, model));
%     out.mcl.sensitivity(model) = mean(results.cl.sensitivity(:, model));
%     out.mcl.Cohenkappa(model) = mean(results.cl.Cohenkappa(:, model));
%     out.mcl.entropy(model) = mean(results.cl.entropy(:, model));
%     
%     out.sDIC(model) = std(results.DIC(:, model));
%     out.scl.chisq(model) = std(results.cl.chisq(:, model));
%     out.scl.phi(model) = std(results.cl.phi(:, model));
%     out.scl.yuleQ(model) = std(results.cl.yuleQ(:, model));
%     out.scl.specificity(model) = std(results.cl.specificity(:, model));
%     out.scl.sensitivity(model) = std(results.cl.sensitivity(:, model));
%     out.scl.Cohenkappa(model) = std(results.cl.Cohenkappa(:, model));
%     out.scl.entropy(model) = std(results.cl.entropy(:, model));
%     
%     for s = 1:sitesize
%         out.mcpt.KL(s,model) = mean(results.cpt.KL{s}(:,model));
%         out.mcpt.entropy(s,model) = mean(results.cpt.entropy{s}(:,model));
%         out.mcpt.accuracy(s,model) = mean(results.cpt.accuracy{s}(:,model));
%         out.mcpt.betadiff(s,model) = mean(results.cpt.betadiff{s}(:,model));
%         out.mcpt.betabwHPD(s,model) = mean(results.cpt.betabwHPD{s}(:,model));
%         out.mcpt.betabwCI(s,model) = mean(results.cpt.betabwCI{s}(:,model));
%         
%         out.scpt.KL(s,model) = std(results.cpt.KL{s}(:,model));
%         out.scpt.entropy(s,model) = std(results.cpt.entropy{s}(:,model));
%         out.scpt.accuracy(s,model) = std(results.cpt.accuracy{s}(:,model));
%         out.scpt.betadiff(s,model) = std(results.cpt.betadiff{s}(:,model));
%         out.scpt.betabwHPD(s,model) = std(results.cpt.betabwHPD{s}(:,model));
%         out.scpt.betabwCI(s,model) = std(results.cpt.betabwCI{s}(:,model));
%     end
end

% fprintf('psi\n')
'psi'
display(mean(results.cl.phi, 1))
display(std(results.cl.phi, [], 1))

'yuleQ'
display(mean(results.cl.yuleQ, 1))
display(std(results.cl.yuleQ, [], 1))

'kappa'
display(mean(results.cl.Cohenkappa, 1))
display(std(results.cl.Cohenkappa, [], 1))

's1'
display(mean(results.cl.sensitivity, 1))
display(std(results.cl.sensitivity, [], 1))

's2'
display(mean(results.cl.specificity, 1))
display(std(results.cl.specificity, [], 1))

'RE'
display(mean(results.cl.relentropy, 1))
display(std(results.cl.relentropy, [], 1))

'SS'
display(results.cl.SS)

'SS2'
display(results.cl.SS2)

mat = [];
for r = 1:REP
    a = [];
    for i = 1:sitesize
        a = [a;results.cpt.accuracy{i}(r,:)];
    end
    mat = [mat; mean(a,1)];
end
'cpt accuracy'
display(mean(mat, 1))
display(std(mat, 1))
% a = reshape(cell2mat(results.cpt.accuracy),[3,49])';
% display(mean(a, 1))

mat = [];
for r = 1:REP
    a = [];
    for i = 1:sitesize
        a = [a;results.fit.DIC{i}(r,:)];
    end
    mat = [mat; mean(a,1)];
end
'fit DIC'
display(mean(mat, 1))
display(std(mat, [], 1))
% a = reshape(cell2mat(results.fit.DIC),[3,49])';
% display(sum(a, 1))

% outH.SS

save('modelCompResults.mat', 'results')
        