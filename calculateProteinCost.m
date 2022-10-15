function result = calculateProteinCost(ecModel,list)
% DESCRIPTION:
% Calculate protein costs for individual reactions in ecModel
% Note that the caluclated protein cost is per flux (mol/gCDW/h)
%
% INPUTS:
% model: an ecModel structure in GECKO2 version
%
% OPTIONAL INPUTS:
% list: a list of reactions with ids consistent with ecModel.rxns
%       (default - all enzymatic reactions in the model)
%
% OUTPUTS:
% result: structure with:
%         * rxns - a list of reactions
%         * costs - calculated protein costs of the corresponding reactions
%                   in the unit gProtein/gCDW per flux (mol/gCDW/h)
%                   NaN means protein cost cannot be calculate due to it
%                   absence or missing kcat and mw data
%
% REFERENCES:
% GECKO2:
% https://doi.org/10.1038/s41467-022-31421-1
% The concept of protein cost
% https://doi.org/10.1073/pnas.1906569116
% https://doi.org/10.1073/pnas.2114622119
%
% Yu Chen, 2022-10-15

% find all enzymes with collected kcat and mw
enzymesIdx = startsWith(ecModel.mets,'prot_') & ~ismember(ecModel.mets,'prot_pool');
enzymesId = ecModel.mets(enzymesIdx);
% identify all enzymatic reactions with collected kcat and mw
listAll = ecModel.rxns(any(ecModel.S(enzymesIdx,:)<0)');

if nargin < 2 || isempty(list)
    list = listAll;
end

result.rxns = list;
result.costs = nan(length(result.rxns),1);

for i = 1:length(list)
    disp([num2str(i) '/' num2str(length(list))]);
    if ismember(list(i),listAll)
        rxnIdx_tmp = ismember(ecModel.rxns,list(i));
        enzymes_tmp = ecModel.mets(ecModel.S(:,rxnIdx_tmp) < 0 & ismember(ecModel.mets,enzymesId));
        [~, b1] = ismember(enzymes_tmp,ecModel.mets);
        coeffs_tmp = -full(ecModel.S(b1,rxnIdx_tmp));
        [~, b2] = ismember(cellfun(@(x) x(6:end),enzymes_tmp,'UniformOutput',false),ecModel.enzymes);
        mws_tmp = ecModel.MWs(b2);
        result.costs(i) = sum(coeffs_tmp.*mws_tmp)*1000; % protein cost = mw/kcat
    else
        warning(['The reaction ' list{i} ' is not in the model or has no enzyme constraint']);
    end
end

