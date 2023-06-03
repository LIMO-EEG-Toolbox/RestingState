function varargout = inter_session_reliability(sessionA,sessionB,metric)

% computes the between reliability based on absolute values and on correlations between session
% for each subject, compute the difference and the Spearman correlation between sessions
%                   compute the differences and the correlations with all other subjects on the other session
% obtain the mean and CI of differences and of correlations
% obtain the frequency of subjects who have a smaller difference or higher
%                         correlation between sessions than among subjects
%
% FORMAT [WithinDiff, HDIW, BetweenDiff, HDIB] = inter_session_reliability(sessionA,sessionB,'metric')
%        [~,~,~,~,freq, indices, opposite, chi, pvalue] = inter_session_reliability(sessionA,sessionB,'metric')
%        
% INPUTS sessionA and sessionB are csv files to be read as table, no raw name expected
%                              measures are computed based on columns
%        'metric' is 'Difference' or 'Spearman' or 'Pearson'
%
% OUTPUTS WithinDiff the trimmed mean difference or correlation of each subject with itself
%         HDIW the highest densiy interval across subjects for Within metric
%         BetweenDiff the trimmed mean difference or correlation of each subject with all other sibjects
%         HDIB the highest densiy interval across subjects for Between metric
%         freq the frequency of subject having a difference smaller than the group or a correlation higher than the group
%         indices flag corrresponding subjects
%         opposite the frequency of between subject metric smaller/bigger than the within subject metric
%         chi is the chi square value evaluated for 1 degree of freedom
%         p-value returns the prob that the observed frequencies are different from 50% 
%         p-value = 1 - chi^2 cumulative distribution function (chi value,1)
%
% sessionA = ['..' filesep 'connectivity_results' filesep 'roi_connect' filesep 'ROI_connect_CS_session-1.csv'];
% sessionB = ['..' filesep 'connectivity_results' filesep 'roi_connect' filesep 'ROI_connect_CS_session-2.csv'];
% [WithinDiff, HDIW, BetweenDiff, HDIB, freq, ~, opposite, chi, pvalue] = inter_session_reliability(sessionA,sessionB,'difference')


%% get the data
if ischar(sessionA)
    sessionA = readtable(sessionA,'delimiter',',');
end

if ischar(sessionB)
    sessionB = readtable(sessionB,'delimiter',',');
end

if any(size(sessionA) ~= size(sessionB))
    error('csv file inputs have different dimensions')
end

%% compute based on differences
if strcmpi(metric,'difference')
    for subject = 1:size(sessionA,2)
        D                    = abs(sessionA{:,subject} - sessionB{:,:});
        if sum(isnan(D(:))) == numel(D)
            WithinDiff(subject) = NaN;
            BetweenDiff(subject) = NaN;
            indices(subject) = 0;
        else
            WithinDiff(subject)  = trimmean(D(:,1),40,'round');
            tmp                  = trimmean(D(:,2:end),40,'round',1);
            BetweenDiff(subject) = trimmean(tmp,40,'round');
            HDI                  = bootCI(tmp);
            indices(subject)     = WithinDiff(subject) < HDI(1) ; % difference smaller than between subjects
        end
    end
    varargout{1} = WithinDiff;
    varargout{2} = bootCI(WithinDiff);
    varargout{3} = BetweenDiff;
    varargout{4} = bootCI(BetweenDiff);
    varargout{5} = mean(indices);
    varargout{6} = indices;
else
    if strcmpi(metric,'Spearman')
        D        = corr(sessionA{:,:},sessionB{:,:},'type','Spearman');
    elseif strcmpi(metric,'Pearson')
        D        = corr(sessionA{:,:},sessionB{:,:},'type','Pearson');
    else
        eror('Only Diffrence, Spearman or Pearson expected as metric')        
    end
    WithinDiff  = diag(D);
    for subject = 1:size(sessionA,2)
        tmp = D(subject,:);
        tmp(subject) = [];
        if sum(isnan(tmp)) == numel(tmp)
            WithinDiff(subject) = NaN;
            BetweenDiff(subject) = NaN;
            indices(subject) = 0;
        else
            BetweenDiff(subject) = trimmean(tmp,40,'round');
            HDI                  = bootCI(tmp);
            indices(subject)     = abs(WithinDiff(subject)) > HDI(2) ; % correlation bigger than between subjects
        end
    end
    varargout{1} = WithinDiff;
    varargout{2} = bootCI(WithinDiff);
    varargout{3} = BetweenDiff;
    varargout{4} = bootCI(BetweenDiff);
    varargout{5} = mean(indices);
    varargout{6} = indices;
end

if nargout > 6
    varargout{7} = mean(varargout{3} < varargout{2}(1));
    opposite     = sum(varargout{3} < varargout{2}(1));
    observed     = [sum(indices) opposite];
    expected     = repmat(size(D,2)/2,1,2);
    varargout{8} = sum((observed-expected).^2./expected);
    varargout{9} = 1-chi2cdf(varargout{7},1,'upper');
end
end

%% Bayes bootstrap
function HDI = bootCI(Data)
prob_coverage = 95/100;
% sample with replacement from Dirichlet
% sampling = number of observations
tmp = sort(Data(~isnan(Data)));
n   = max(size(tmp)); bb = zeros(1000,1);
for boot=1:1000 % bootstrap loop
    theta    = exprnd(1,[n,1]);
    weigths  = theta ./ repmat(sum(theta),n,1);
    resample = datasample(tmp,n,'Replace',true,'Weights',weigths);
    bb(boot) = trimmean(resample,40,'round');
end
sorted_data   = sort(bb); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,1)); % upper bound
nCIs          = size(sorted_data,1) - upper_centile;
ci            = 1:nCIs;
ciWidth       = sorted_data(ci+upper_centile) - sorted_data(ci); % all centile distances
[~,index]     = find(ciWidth == min(ciWidth)); % densest centile
if length(index) > 1 % many similar values
    index = index(1);
end
HDI          = [sorted_data(index) sorted_data(index+upper_centile)];
end


