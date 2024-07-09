clc;

% Define the cell array
C = { [1 2 3] [3 6 9]};

% Calculate the average of all cell elements
avgElement = cellArrayAverage(C);
disp('The average of all cell elements is:');
disp(avgElement);

function avg = cellArrayAverage(C)
    % Initialize the sum array with zeros of the same size as the first element
    sumArray = zeros(size(C{1}));
    
    % Initialize counter
    numElements = numel(C);
    
    % Loop through each cell and sum the elements
    for i = 1:numElements
        if isnumeric(C{i}) && isequal(size(C{i}), size(sumArray))
            sumArray = sumArray + C{i};
        else
            error('All elements must be numeric arrays of the same size.');
        end
    end
    
    % Calculate the average
    avg = sumArray / numElements;
end