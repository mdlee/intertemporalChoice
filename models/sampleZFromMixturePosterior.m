function z = sampleZFromMixturePosterior(Z, nModels)
%SAMPLEZFROMMIXTUREPOSTERIOR  One z vector; each person from their marginal.

nP = size(Z, 2);
nDraws = size(Z, 1);
z = zeros(nP, 1);
for i = 1:nP
  z(i) = Z(randi(nDraws), i);
end
z = max(1, min(nModels, round(z(:))));
end
