function inter_session_reliability_test

% tests that inter_session_reliability.m outputs what it is supposed too
% under the null, all data are randn generated
% under perfect condition, all data are well correlated, with 50% subject
%                          more reliables 

[WithinDiff, HDIW, BetweenDiff, HDIB, freq, indices] = inter_session_reliability(array2table(randn(100,1000)),...
    array2table(randn(100,1000)),'difference');
if freq ~= 0
    error('0% of within subject expecxted diffeent')
end

[WithinDiff, HDIW, BetweenDiff, HDIB, freq, indices] = inter_session_reliability(array2table(randn(100,1000)),...
    array2table(randn(100,1000)),'Spearman');
