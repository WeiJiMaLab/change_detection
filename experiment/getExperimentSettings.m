% This is a helper function to set the settings for the experiment code. It
% takes in the subject initials, and could be modified to generate
% experimental trials in a pseudorandom (rather than random according to
% the generative model) way.

function settings = getExperimentSettings()

% some initializations
exptype='';

deltatype='';
deltaval = [];
deltavartype = '';
deltavartypestd=[];
deltavartypemu=[];
deltavartyperange = [];

reliabilitytype='';
reliabilityval=[];
reliabilityvartype='';
reliabilityvartypestd=[];
reliabilityvartypemu=[];

setsizeval=[];

multiplicity = [];

drawstimuli = '';

nTrials = [];

degreeres = [];

variableparam = '';

breaknum =[];

expID = '';
clc;

% Get subject ID
subjid   = input('Subject initials (three letters)......: ','s');

% Only do change detection experiment w/ different set sizes
expID = 'CDSetSize';

if(~strcmp(expID,'CDSetSize') && ~strcmp(expID,'Training'))
    % Threshold or actual experiment? How many trials per condition?
    while ~(strcmp(exptype,'threshold') || strcmp(exptype,'actual'))
        exptype = input('Experiment type (threshold/actual)....................: ','s');
    end
    while ~(strcmp(drawstimuli,'uniformly') || strcmp(drawstimuli,'randomly'))
        drawstimuli = input('How to draw stimuli (uniformly/randomly)..............: ','s');
    end
    while ~(strcmp(exptype,'threshold') || ~isempty(multiplicity) || strcmp(drawstimuli,'randomly'))
        multiplicity = input('Trials per condition (integer)........................: ','s');
    end
    while ~(~isempty(nTrials) || strcmp(drawstimuli,'uniformly'))
        nTrials = input('How many trials (integer)....................: ','s');
    end
    while ~(~isempty(degreeres))
        degreeres = input('How fine of angle sampling (degrees)..................: ','s');
    end
    
    % What type of delta? What parameters for the type?
    while ~(strcmp(deltatype,'constants') || strcmp(deltatype,'variable'))
        deltatype = input('Delta type (constants/variable).......................: ','s');
    end
    while ~(strcmp(deltatype,'variable') || ~isempty(deltaval))
        deltaval = input('Delta values (vector double [a b ...])................: ','s');
    end
    while ~(strcmp(deltatype,'constants') || strcmp(deltavartype,'Uniform') || strcmp(deltavartype,'Gaussian') || strcmp(deltavartype,'Von_Mises'))
        deltavartype = input('Delta distribution (Gaussian/Von_Mises/Uniform).......: ','s');
    end
    while ~(strcmp(deltavartype,'Uniform') || ~isempty(deltavartyperange))
        deltavartyperange = input('Delta range (0:degrees from 0-180).....................: ','s');
    end
    while ~(strcmp(deltavartype,'Uniform') || strcmp(deltatype,'constants') || ~isempty(deltavartypestd))
        deltavartypestd = input('Delta standard deviations (vector double [a b ...])...: ','s');
    end
    while ~(strcmp(deltavartype,'Uniform') || strcmp(deltatype,'constants') || ~isempty(deltavartypemu))
        deltavartypemu = input('Delta means (mu) (vector double [a b ...])............: ','s');
        s1 = size(deltavartypestd);
        s2 = size(deltavartypemu);
        if(s1(1) ~= s2(1) || s1(2) ~= s2(2))
            fprintf('Vectors of standard deviation and mean must be the same size.\n')
            deltavartypemu=[];
        end
    end
    
    % What type of reliability? What parameters?
    while ~(strcmp(reliabilitytype,'constants') || strcmp(reliabilitytype,'variable'))
        reliabilitytype = input('Reliability type (constants/variable).................: ','s');
    end
    while ~(strcmp(reliabilitytype,'variable') || ~isempty(reliabilityval))
        reliabilityval = input('Reliability value (vector double [a b ...])...........: ','s');
    end
    while ~(strcmp(reliabilitytype,'constants') || strcmp(reliabilityvartype,'Gaussian') || strcmp(reliabilityvartype,'Von_Mises'))
        reliabilityvartype = input('Reliability distribution (Gaussian/Von_Mises).........: ','s');
    end
    while ~(strcmp(reliabilitytype,'constants') || ~isempty(reliabilityvartypestd))
        reliabilityvartypestd = input('Reliability standard deviations (vector double [a b ...]): ','s');
    end
    while ~(strcmp(reliabilitytype,'constants') || ~isempty(reliabilityvartypemu))
        reliabilityvartypemu = input('Reliability means (mu) (vector double [a b ...])......: ','s');
        s1 = size(reliabilityvartypestd);
        s2 = size(reliabilityvartypemu);
        if(s1(1) ~= s2(1) || s1(2) ~= s2(2))
            fprintf('Vectors of standard deviation and mean must be the same size.\n')
            reliabilityvartypemu=[];
        end
    end
    
    % What set sizes?
    while ~(~isempty(setsizeval))
        setsizeval = input('Set sizes (vector integer [a b ...])..................: ','s');
    end
    
    % Which parameter to vary within block
    while ~(strcmp(variableparam,'delta') || strcmp(variableparam,'setsize') || strcmp(variableparam,'reliability'))
        variableparam = input('Parameter to vary in block (delta/setsize/reliability): ','s');
    end
    
    % How many breaks?
    while ~(~isempty(breaknum))
        breaknum = input('Number of breaks: ','s');
    end
elseif(strcmp(expID,'CDSetSize'))
    % Set parameters for the first change detection experiment, variable
    % set size per block
    exptype='actual';
    
    deltatype='variable';
    deltavartyperange = 180;
    deltavartype = 'Uniform';
    
    reliabilitytype='constants';
    reliabilityval= .9;
    
    setsizeval=[2 4 6 8];    
    multiplicity = 75;
    drawstimuli = 'randomly';   
    nTrials = 600;
    degreeres = 1;
    variableparam = '';
    breaknum = 6;
    
    
elseif(strcmp(expID,'Training'))
    % Set parameters for the first change detection experiment, variable
    % set size per block
    exptype='train';
    
    deltatype='variable';
    deltavartyperange = 180;
    deltavartype = 'Uniform';
    
    reliabilitytype='constants';
    reliabilityval= .9;
    
    setsizeval=[2 4 6 8];    
    multiplicity = 75;
    drawstimuli = 'randomly';   
    nTrials = 16;
    degreeres = 1;
    variableparam = '';
    breaknum = 4;
end

settings.subjid = upper(subjid);
settings.exptype = exptype;

settings.deltatype = deltatype;
settings.deltaval = deltaval;
settings.deltavartype = deltavartype;
settings.deltavartypestd = deltavartypestd;
settings.deltavartypemu = deltavartypemu;
settings.deltavartyperange = deltavartyperange;

settings.reliabilitytype = reliabilitytype;
settings.reliabilityval = reliabilityval;
settings.reliabilityvartype = reliabilityvartype;
settings.reliabilityvartypestd = reliabilityvartypestd;
settings.reliabilityvartypemu = reliabilityvartypemu;

settings.expID = expID;

settings.setsizeval = setsizeval;

settings.multiplicity = multiplicity;

settings.degreeres = degreeres;

settings.variableparam = variableparam;

settings.breaknum = breaknum;

settings.drawstimuli = drawstimuli;

settings.stimtime = 100;

if(strcmp(drawstimuli,'uniformly'))
    settings.nTrials = 2*length(deltaval)*length(reliabilityval)*length(setsizeval)*multiplicity;
elseif(strcmp(drawstimuli,'randomly'))
    settings.nTrials = nTrials;
end
