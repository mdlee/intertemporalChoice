function generateHierarchicalMuPrecJags(srcJagsPath, dstJagsPath, precScale)
%GENERATEHIERARCHICALMUPRECJAGS  Scale dnorm precisions on group mu_* hyperpriors.
%
%   PRECSCALE = 0.5 halves prior precision (wider priors on group means);
%   PRECSCALE = 2 doubles prior precision (tighter priors on group means).
%
%   Works for latent-mixture and separate hierarchical cognitive models.

txt = fileread(srcJagsPath);
lines = splitlines(string(txt));
if isempty(lines)
  error('Empty JAGS file: %s', srcJagsPath);
end

nChanged = 0;
for i = 1:numel(lines)
  line = char(lines(i));
  tok = regexp(line, '^\s*(mu_\w+)\s*~\s*dnorm\(([^,]+),\s*([^)]+)\)(.*)$', 'tokens', 'once');
  if isempty(tok)
    continue;
  end
  prec = str2double(strtrim(tok{3}));
  if ~isfinite(prec)
    error('Could not parse precision on line %d: %s', i, line);
  end
  newPrec = prec * precScale;
  lines(i) = sprintf('  %s ~ dnorm(%s, %s)%s', tok{1}, strtrim(tok{2}), ...
    formatJagsPrec(newPrec), tok{4});
  nChanged = nChanged + 1;
end

if nChanged == 0
  error('generateHierarchicalMuPrecJags:noMuPriors', ...
    'No mu_* ~ dnorm(...) lines found in %s', srcJagsPath);
end

header = sprintf('# Generated from %s; mu_* dnorm precision × %g\n', ...
  srcJagsPath, precScale);
outTxt = char(strjoin([string(header); lines], newline));
fid = fopen(dstJagsPath, 'w');
if fid < 0
  error('Could not write %s', dstJagsPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', outTxt);
fprintf('Wrote %s (mu prior precision × %g; %d lines)\n', dstJagsPath, precScale, nChanged);
end

function s = formatJagsPrec(p)
if abs(p - round(p)) < 1e-9 && p >= 1 && p < 1e6
  s = sprintf('%d', round(p));
elseif p >= 1e-3
  s = sprintf('%.4g', p);
else
  s = upper(sprintf('%.1E', p));
end
end
