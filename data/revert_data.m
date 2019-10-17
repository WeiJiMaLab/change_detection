function subjCell_new = revert_data(subjCell)

% Code for reverting data from the intuitive format to what the code uses

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
    
    data = subjCell{ii};
    
    curr_data = zeros(length(data.change),71);
    curr_data(:,1) = double(data.change);
    curr_data(:,2) = double(data.response);
    curr_data(:,4) = double(data.RT);
    curr_data(:,5) = double(data.setsize);
    curr_data(:,56:63) = double(round((180/pi)*data.ort_first));
    curr_data(:,64:71) = double(round((180/pi)*data.ort_second));
    
    for jj = 1:size(curr_data,1)
        curr_data(jj,(56+data.setsize(jj)):63)=1;
        curr_data(jj,(64+data.setsize(jj)):71)=1;
    end
    
    subjCell_new{ii} = curr_data;
    
end
        