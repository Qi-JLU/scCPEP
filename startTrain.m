% take multiclass classification and voting for combining base classifiers as an example.
function results = startTrain()

data = load([pwd, filesep, 'ensemble.mat']);

global PredictClass TrueClass

PredictClass = double(data.PredictClass);
TrueClass = double((data.TrueClass)');


% The objective function f is taken as the error on the validation data set. 
% solution: a subset of base classifiers, i.e., a Boolean string of n; 
function fValue=f(solution)
    subPredInd = find(solution == 1);
    subLength = length(subPredInd);
    subMat = PredictClass(:, subPredInd);
    tag = 1./TrueClass;
    res = subMat .* tag; % if pred==true, the value will be 1
    result = sum(res==1, 2); % row
    
    if subLength>0
        fValue= sum((2*result)<subLength)+sum((2*result)==subLength)/2;
    else
        fValue=inf;
    end
end

% The eval function is also taken as the error on the validation data set.
function evalValue=eval(solution)    
    subPredInd = find(solution == 1);
    subLength = length(subPredInd);
    subMat = PredictClass(:, subPredInd);
    tag = 1./TrueClass;
    res = subMat .* tag; % if pred==true, the value will be 1
    result = sum(res==1, 2); % row
    
    if subLength>0
        evalValue=sum((2*result)<subLength)+sum((2*result)==subLength)/2;
    else
        evalValue=inf;
    end
end

n = size(PredictClass, 2);

% use the PEP function to select the subEnsemble.
% results = PEP(n,@f,@eval);
results = PEP_acceleration(n,@f,@eval);

end
