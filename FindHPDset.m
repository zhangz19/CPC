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