% parse intertemporal choice data

clear;

d.info = {'Elena data'};
d.nParticipants = 25;
d.nTrials = 240;
d.nPairs = 24;
for i = 1:d.nParticipants
  T = readtable(sprintf('sujeto_%d_tiempo.csv', i), 'delimiter', ',');
  for j = 1:d.nTrials
    d.trial(i, j) = T{j, 1};
    d.rewardA(i, j) = T{j, 2};
    d.timeA(i, j) = T{j, 3};
    d.rewardB(i, j) = T{j, 4};
    d.timeB(i, j) = T{j, 5};
    [~, d.decision(i, j)] = ismember(T{j, 6}, {'left', 'right'});
    d.pair(i, j) = T{j, 8};
    d.rt(i, j) = T{j, 9};
    d.LLchosen(i, j) = T{j, 7};
    if d.rewardA(i, j) > d.rewardB(i, j)
      d.LL(i, j) = (d.decision(i, j) == 1);
    else
      d.LL(i, j) = (d.decision(i, j) == 2);
    end
  end
end

% order participants
[~, d.participantOrder] = sort(mean(d.LL, 2));
d.trialOrder = nan(d.nParticipants, d.nTrials);
d.trial = d.trial(d.participantOrder, :);
d.rewardA = d.rewardA(d.participantOrder, :);
d.rewardB = d.rewardB(d.participantOrder, :);
d.timeA = d.timeA(d.participantOrder, :);
d.timeB = d.timeB(d.participantOrder, :);
d.decision = d.decision(d.participantOrder, :);
d.pair = d.pair(d.participantOrder, :);
d.rt = d.rt(d.participantOrder, :);
d.LL = d.LL(d.participantOrder, :);
d.LLchosen = d.LLchosen(d.participantOrder, :);

% classify pairs
d.pairTypes = {{'small', 'intervals'}, {'medium', 'intervals'}, ...
   {'large', 'intervals'}, {'small', 'outcomes'}};
d.pairType = nan(d.nParticipants, d.nTrials);
d.nPairTypes = numel(d.pairTypes);
for i = 1:d.nParticipants
   for j = 1:d.nTrials
      if d.rewardA(i, j) < 2000 % small outcomes are all 1000-2000 pesos
         d.pairType(i, j) = 4;
      else
         if abs(d.timeA(i, j)-d.timeB(i, j)) == 1
            d.pairType(i, j) = 1;
         elseif abs(d.timeA(i, j)-d.timeB(i, j)) == 2
            d.pairType(i, j) = 2;
         else
            %abs(d.timeA(i, j)-d.timeB(i, j))
            d.pairType(i, j) = 3;
         end
      end
   end
end

d.LLtrial = nan(size(d.LL));
for i = 1:d.nParticipants
   tmp = [];
   for k = 1:d.nPairTypes
      tmp = [tmp find(d.pairType(i, :) == k)];
   end
   d.trialOrder(i, :) = tmp;
   d.LLtrial(i, :) = d.LL(i, tmp);
   d.LLchosen(i, :) = d.LLchosen(i, tmp);
end
d.LL = double(d.LL);

save intertemporalChoice d

