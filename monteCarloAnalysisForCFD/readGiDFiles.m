function fileList = readGiDFiles()
    % fileList = readGiDFiles()
    % Return a cell array containing every GiD .dat file from
    % root/inputGiD/FEMComputationalFluidDynamicsAnalysis

    % Read folder contents
    fileNames = dir('../inputGiD/FEMComputationalFluidDynamicsAnalysis');
    fileNames = {fileNames.name};
    % Check file names
    fileList = cell.empty(0,1);
    for k=1:length(fileNames)
        if length(fileNames{k})>3 && strcmp(fileNames{k}(end-3:end), '.dat')
            fileList{end+1} = fileNames{k}(1:end-4);
        end
    end

end