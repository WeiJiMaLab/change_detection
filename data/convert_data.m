function subjCell_new = convert_data(subjCell)

% Code for converting data into a more intuitive format

% Reference
% 1-> target presence (1,0)
% 2-> response (1,0)
% 3-> real stim time (ms)
% 4-> reaction time (ms)
% 5-> number of stimuli [setsize] (int)
% 6-> delta
% 56:63-> first array orientations (radians)
% 64:71-> second array orientations (radians)

subjCell_new = cell(size(subjCell));

for ii = 1:length(subjCell)
        
    curr_data = subjCell{ii};
    
    data.change = curr_data(:,1)>0;
    data.response = curr_data(:,2)>0;
    temp_delta = .5*sum(abs(circ_dist((pi/90)*curr_data(:,56:63),(pi/90)*curr_data(:,64:71))),2);
    temp_delta(temp_delta < eps) = 0; % numerical error in computing difference
    data.delta = temp_delta;
    data.RT = curr_data(:,4);
    data.setsize = curr_data(:,5);
    data.ort_first = curr_data(:,56:63);
    data.ort_second = curr_data(:,64:71);
    
    % convert to radians, making sure to leave out the non-stimuli that
    % were set to 1 for the experiment
    for jj = 1:size(curr_data,1)
        data.ort_first(jj,(curr_data(jj,5)+1):end)=0;
        data.ort_second(jj,(curr_data(jj,5)+1):end)=0;
    end
    data.ort_first = data.ort_first*(pi/180);
    data.ort_second = data.ort_second*(pi/180);
    
    subjCell_new{ii} = data;
    clear data
    
end
        