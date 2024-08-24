% Open the original file
originalFile = fopen('Fe2O3Pt_04_40to142K.txt', 'r');

% Check if the file was opened successfully
if originalFile == -1
    error('Could not open the original file.');
end

try
    x = 24831;
    % Initialize the line counter and file counter
    lineCounter = 0;
    fileCounter = 1;

    % Open the first output file
    outputFile = fopen(sprintf('Fe2O3Pt_04_%d.txt', fileCounter), 'w');

    % Check if the file was opened successfully
    if outputFile == -1
        error('Could not open the output file.');
    end

    % Read the original file line by line
    while ~feof(originalFile)
        line = fgetl(originalFile);
        fprintf(outputFile, '%s\n', line);
        lineCounter = lineCounter + 1;

        % If we have written x lines to the current output file
        if mod(lineCounter, x) == 0
            % Close the current output file
            fclose(outputFile);

            % Increment the file counter
            fileCounter = fileCounter + 1;

            % Open the next output file
            outputFile = fopen(sprintf('Fe2O3Pt_04_%d.txt', fileCounter), 'w');

            % Check if the file was opened successfully
            if outputFile == -1
                error('Could not open the output file.');
            end
        end
    end

    % Close the last output file
    fclose(outputFile);
catch ME
    % If an error occurs, close any open files before rethrowing the error
    fclose('all');
    rethrow(ME);
end

% Close the original file
fclose(originalFile);
