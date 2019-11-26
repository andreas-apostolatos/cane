function inputList = readInputFunctions()
    % inputList = readInputFunctions()
    % Return a cell array containing every function name that matches
    % input_* from root/inputFunctions

    % Read folder contents
    fileNames = dir('inputFunctions');
    fileNames = {fileNames.name};
    % Check file names
    inputList = cell.empty(0,1);
    for k=1:length(fileNames)
        if length(fileNames{k})>6
            if strcmp( fileNames{k}(1:6), 'input_' )
                inputList{end+1} = fileNames{k}(1:end-2);                   % <--- Chop .m from the end
            end
        end
    end

end