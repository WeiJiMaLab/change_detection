function [P_C_HAT_Cell] = load_ML_fits(model_num,is_color)

% choose which type of data to load. 0 and 1 are the relevant ones
switch is_color
    case 0
        load ../data/subjCell
        num_subjects = length(subjCell);
        prefix = [];
    case 1
        load ../data/subjColorCell
        num_subjects = length(subjColorCell);
        prefix = 'Color_';
    case 2
        load ./subjFakeCell
        num_subjects = length(subjFakeCell);
        prefix = 'Fake_';
    case 3
        load ./csubjFakeCell
        num_subjects = length(subjFakeCell);
        prefix = 'cFake_';
    case 4
        load ./sasubjFakeCell
        num_subjects = length(subjFakeCell);
        prefix = 'saFake_';
    case 5
        load ./csasubjFakeCell
        num_subjects = length(subjFakeCell);
        prefix = 'csaFake_';
        
end

    P_C_HAT_Cell = cell(num_subjects,1);
    
    for i = 1:num_subjects
        try
            load(['max_P_C_HAT/' prefix 'max_P_C_HAT_' num2str(i) '_' num2str(model_num)])
            P_C_HAT_Cell{i} = max_P_C_HAT;
        catch
            fprintf('\nCould not load subject %g!\n',i)
        end
    end

        
end
