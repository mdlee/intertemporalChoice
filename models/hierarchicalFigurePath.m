function pathOut = hierarchicalFigurePath(figuresDir, modelName, tag, doSave)
%HIERARCHICALFIGUREPATH  PNG path under models/figures for hierarchical fits.
if doSave
  if ~isfolder(figuresDir)
    mkdir(figuresDir);
  end
  pathOut = fullfile(figuresDir, sprintf('%s_%s.png', modelName, tag));
else
  pathOut = '';
end
end
