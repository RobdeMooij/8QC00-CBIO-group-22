function sysPack = getDualSysPack(systemTemplate,modVars,randmodulator)
% find all combinations
fn = fieldnames(modVars);

% create parameter sets for FIRST RUN
myCount = 1;
for ci=1:numel(modVars.(fn{1}))
    for cj=1:numel(modVars.(fn{2}))
         % write system variables to current set, vary random seeds
          sysPack{ci,cj} = systemTemplate;
          sysPack{ci,cj}.params.initialSeed = myCount * randmodulator;
          myCount = myCount+1;
    end
end

sysPack = modifyTwoParams(sysPack,modVars);

end