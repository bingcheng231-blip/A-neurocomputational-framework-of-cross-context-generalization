function [negLogLik,VV,PP, step_lr] = s003_ctx_RL_LogLikelihood(theta,stimuli_id,outcomes,choices,STRUCT_flag, decayed_lr_flag, prior_structure)
%   s003_ctx_RL_LogLikelihood 此处显示有关此函数的摘要
%   此处显示详细说明

    nStimuli = length(unique(stimuli_id)); 

    % learning rate setting
    lr= theta(1); % initial learning rate
    decayed_lr = lr;

    if decayed_lr_flag
        decay_rate = theta(2);
        decay_steps = theta(3);
        beta = theta(4);
    else
        beta = theta(2);
    end

    if STRUCT_flag
        cross_term = squareform(theta(end-275:end));
    else % NAIVE model
        cross_term = zeros(nStimuli, nStimuli);
    end

    if ~isempty(prior_structure)
        if ~STRUCT_flag
            middle_point = theta(end);
            prior_structure = prior_structure-middle_point;
            max_ct = max(max(prior_structure), abs(min(prior_structure)));
            prior_structure = prior_structure./max_ct;
            cross_term = squareform(prior_structure);
        else
            error('Structure estimation is in conflict with the prior structure')
        end
    end

    v0 = 0; % initial value at beginning of block
    v = nan(nStimuli,1);
    v(:)=v0;

    Trials_num = length(stimuli_id); % number of trials
    PP = nan(Trials_num,1); % track prob of accepting through time 
                    % (pass value after last trial through sigmoid).  this is
                    % between 0 and 1). This is used before viewing the current
                    % outcome. 
    VV = nan(Trials_num, nStimuli); % track value through time,this is between -1 and 1.
                                    % This is the updated value, after viewing the outcome.
    step_lr = nan(Trials_num,1);

    for ti = 1:Trials_num % trials   
        sid = stimuli_id(ti); % current stimulus   
        % v is the EV table
        p = 1/(1+exp(- beta*(v(sid)))); % pass through sigmoid to obtain probabilities between 0 and 1. 
        PP(ti) = p; % probability for GOOD
        
        % update value according to the current trial outcome
        dv = outcomes(ti) - v(sid);
        v(sid) = v(sid) + decayed_lr*dv;
        
        % update other value according to their correlation to the current trial outcome
        % This part of the updating is available only for the structural model.
        other_ccn = setdiff(1:nStimuli, sid);
        ct = cross_term(other_ccn, sid); % cross term is an object relation matrix (stimulus num x stimulus num)
        v(other_ccn) = (1 - abs(ct).*decayed_lr) .* v(other_ccn) + ct.*decayed_lr.*outcomes(ti);  
        VV(ti,:) = v;   

        if decayed_lr_flag
            step_lr(ti) = decayed_lr;
            decayed_lr = lr.*(decay_rate.^(ti./decay_steps));
        end
    end
    PP(PP==0) = 1e-16;
    PP(PP==1) = 1-1e-16;
    negLogLik = sum(-log(1-PP(choices==0)))+sum(-log(PP(choices==1)));
end

















































