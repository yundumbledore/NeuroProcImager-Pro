function main(data_file, tasks)
    if ismember('parameter estimation', tasks)
        disp(['Now running Parameter Estimation on recording #' num2str(data_file) '...'])

        addpath(genpath('./inference'))
        inference(data_file)
        rmpath(genpath('./inference'))
    end

    if ismember('stability analysis', tasks)
        disp(['Now running Stability Analysis on recording #' num2str(data_file) '...'])

        addpath(genpath('./stability_analysis'))
        % stability_analysis(data_file)
        show_dynamic_stability(data_file)
        rmpath(genpath('./stability_analysis'))
    end

    if ismember('chaos analysis', tasks)
        disp(['Now running Chaos Analysis on recording #' num2str(data_file) '...'])
        addpath(genpath('./chaos_analysis'))
        chaos_analysis(data_file)
        show_dynamic_chaos(data_file)
        rmpath(genpath('./chaos_analysis'))
    end
end