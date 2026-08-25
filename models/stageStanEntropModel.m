function [modelPath, workDir] = stageStanEntropModel(modelsDir, modelName)
%STAGESTANENTROPMODEL  Inline include/entrop.stan for Trinity/CmdStan (no #include).
%
%   Robustness jobs are named ...MuPrecHalf_entrop / MuPrecDouble_entrop
%   but share the baseline .stan; muPrecScale is passed in the data struct.

stanDir = fullfile(modelsDir, 'stan');
baseName = regexprep(char(modelName), 'MuPrec(Half|Double)_', '');
src = fullfile(stanDir, [baseName '.stan']);
if ~isfile(src)
  error('stageStanEntropModel:missingModel', 'Stan model not found: %s', src);
end
workDir = fullfile(stanDir, 'tmp', char(modelName));
if ~isfolder(workDir)
  mkdir(workDir);
end
modelText = fileread(src);
incPath = fullfile(stanDir, 'include', 'entrop.stan');
if contains(modelText, '#include') && isfile(incPath)
  incText = fileread(incPath);
  modelText = regexprep(modelText, '#include\s+"entrop\.stan"\s*', '');
  modelText = sprintf('%s\n%s', incText, modelText);
end
modelPath = fullfile(workDir, [baseName '.stan']);
fid = fopen(modelPath, 'w');
if fid < 0
  error('stageStanEntropModel:write', 'Cannot write %s', modelPath);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s', modelText);
end
