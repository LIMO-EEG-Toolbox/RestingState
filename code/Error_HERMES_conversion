export_folder = fullfile(pwd, ['H_connectivity_results' filesep 'HERMES_CSVS']);
hermes_metrics = {'COH', 'COR', 'iCOH','DPI','PLI', 'PLV','RHO','wPLI'};
mkdir(export_folder);

num_subs = size(indexes.(hermes_metrics{1}).data, 2);
num_sessions = size(indexes.(hermes_metrics{1}).data, 1);

% Check if a given metric needs to be averaged over frequencies/bands by checking and conditioning on each metrics' dimensions:
for s = 1:length(hermes_metrics)
    if (all(~any(strcmp(indexes.(hermes_metrics{s}).dimensions, 'frequency')))  == 0 | all(~any(strcmp(indexes.(hermes_metrics{s}).dimensions, 'band')))  == 0)
        for i = 1:num_subs
            for j = 1:num_sessions
                data = indexes.(hermes_metrics{s}).data{j, i}; % Getting the data from the metric
                averagedData = mean(data, 3); % Compute the average along the third dimension (that is frequencies/bands)
                indexes.(hermes_metrics{s}).data{j, i} = averagedData; % Assign the averaged data back to the metrics such that we have averaged if needed and all metrics have same dimensions
            end
        end
    end
end



%Convertion to csv files
for s = 1:length(hermes_metrics) 
    metric_data = indexes.(hermes_metrics{s}).data;
    for j = 1:num_sessions
        session_data = []; % Pre allocate to store data for each session
        
        for i = 1:num_subs
            % Symmetric data matrix
            data = triu(metric_data{j, i});
            data_vector = data(triu(true(size(data)), 1));
            % Append the data vector for the subject into the corresponding session_data
            session_data = [session_data, data_vector];
        end
        
        table_data = array2table(session_data);
        var_names = cellstr(num2str((1:num_subs)', 'Subject_%d'));
        table_data.Properties.VariableNames = var_names;
        
        % Create csv file for the current session:
        csv_filename = fullfile(export_folder, ['session' num2str(j) '_' hermes_metrics{s} '.csv']);
        writetable(table_data, csv_filename);
    end
end



