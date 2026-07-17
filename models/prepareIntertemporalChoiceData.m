function [data, d] = prepareIntertemporalChoiceData(dataName, dataDir)
%PREPAREINTERTEMPORALCHOICEDATA  Load data and build LL/SS-oriented trial matrices.
if nargin < 2 || isempty(dataDir)
  dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
end
if nargin < 1 || isempty(dataName)
  dataName = 'intertemporalChoice';
end

load(fullfile(dataDir, dataName), 'd');

decision = zeros(d.nParticipants, d.nTrials);
rLL = zeros(d.nParticipants, d.nTrials);
rSS = zeros(d.nParticipants, d.nTrials);
tLL = zeros(d.nParticipants, d.nTrials);
tSS = zeros(d.nParticipants, d.nTrials);

for i = 1:d.nParticipants
  for j = 1:d.nTrials
    if d.rewardA(i, j) > d.rewardB(i, j)
      decision(i, j) = (d.decision(i, j) == 1);
      rLL(i, j) = d.rewardA(i, j);
      rSS(i, j) = d.rewardB(i, j);
      tLL(i, j) = d.timeA(i, j);
      tSS(i, j) = d.timeB(i, j);
    else
      decision(i, j) = (d.decision(i, j) == 2);
      rLL(i, j) = d.rewardB(i, j);
      rSS(i, j) = d.rewardA(i, j);
      tLL(i, j) = d.timeB(i, j);
      tSS(i, j) = d.timeA(i, j);
    end
  end
end

data = struct( ...
  'nParticipants', d.nParticipants, ...
  'nTrials', d.nTrials, ...
  'rLL', rLL, ...
  'rSS', rSS, ...
  'tLL', tLL, ...
  'tSS', tSS, ...
  'decision', decision);

end
