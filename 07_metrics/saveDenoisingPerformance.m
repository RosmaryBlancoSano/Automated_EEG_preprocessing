
function saveDenoisingPerformance(results, fpath_out, channelConfig)
    % Unique participants from the field names of results
    fieldNames = fieldnames(results);
    subjects = unique(cellfun(@(x) strtok(x, '_'), fieldNames, 'UniformOutput', false));

    % One folder per participant
    for i = 1:length(subjects)
        subj = subjects{i};  % Es. "Subj1", "Subj2"

        % Create the folder of the participant if it does not exist
        subj_folder = fullfile(fpath_out, subj);
        if ~exist(subj_folder, 'dir')
            mkdir(subj_folder);
        end

        % Cell array for the data
        dataCell = {};
        header = {'Subject', 'Method', 'ElapsedTime', 'CPUTime', 'MemoryUsed'};

        % Methods of the current participant
        subj_fields = fieldNames(startsWith(fieldNames, subj));

        for j = 1:length(subj_fields)
            fieldname = subj_fields{j};
            
            % Method name (for example "EMD" from "Subj1_EMD")
            method = strrep(fieldname, [subj '_'], '');

            % Performance values
            elapsedTime = results.(fieldname).elapsedTime;
            cpuTime = results.(fieldname).cpuTime;
            memoryUsed = results.(fieldname).memoryUsed;

            % Add the data to the cell array
            dataCell = [dataCell; {subj, method, elapsedTime, cpuTime, memoryUsed}];
        end

        % Convert the cell array to a table
        T = cell2table(dataCell, 'VariableNames', header);

        % Name of the CSV file of the participant
        
        csv_name = sprintf('denoising_performance_%s_%s.csv', subj, channelConfig);
        writetable(T, fullfile(subj_folder, csv_name));

        fprintf('Results saved for Subj%d in %s\n', subj, csv_name);
    end
end

