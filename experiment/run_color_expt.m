% Run the COLOR experiment as run in the paper, generating trials from the
% generative model outlined in the supplementary figures. Note: Make sure
% that the most recent version of Psychtoolbox is installed. This one runs
% well on Psychtoolbox 3.0.10. This will crash if multiple monitors are
% used.

% There will first be an instruction screen, then the trials will commence
% immediately. There are untimed breaks throughout the experiment. The
% subject will use the "y" and "u" keys to respond. To exit the experiment
% prematurely, press "esc". Each run is 600 trials, with feedback on each
% trial. 

% IMPORTANT: Check column 3 of the data file in the output directory for
% the actualy stimulus display time. This can change depending on the
% monitor refresh rate. If it is not 117, adjust the "settings.stimtime"
% property in getExperimentSettings so that it is. It's best to do this as
% a test before running the experiment. We recommend using a monitor with a
% 60 Hz refresh rate.

% This will save an output file in the "output" directory that contains the
% initials of the subject as well as an identifier of the time and date the
% session was run. Each separate run of the experiment generates a
% different output file, which should be appended together into a data
% structure (see analysis README) before analyzing.

function run_color_expt()

% Set random stream
rng('shuffle')

try
    
    % Initialization
    clear all;

    % Get settings from experimenter
    settings = getExperimentSettings();
    multiplicity = settings.multiplicity;
    setsizeval = settings.setsizeval;
    deltatype = settings.deltatype;
    deltavartype = settings.deltavartype;
    deltaval = settings.deltaval;
    reliabilityval = settings.reliabilityval;
    deltavartyperange = settings.deltavartyperange*(pi/180);
    degreeres = settings.degreeres;
    
    if isempty(settings.nTrials)
        nTrials = 2*multiplicity*length(setsizeval)*length(deltaval)*length(reliabilityval);
    else
        nTrials = settings.nTrials;
    end
    
    breaknum = settings.breaknum;
    
    drawstimuli  = settings.drawstimuli;
    
    subjid = settings.subjid;
    settings.stimtype = 'ellipse';
    settings.disttype = 'hetero';
    settings.stimsize = .4;
    
    % set some condition-independent variables
    settings.makeScreenShot  = 0;    % if 1, then screenshots of stimuli will be made.
    % IMPORTANT: If this is set to 1, the stimuli will be shown for WAY TOO
    % LONG (i.e. 1 second or more). ONLY use this to get screenshots for
    % viewing later, not when actually taking data.

    settings.screen_width    = 40;   % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)
    settings.barwidth        = .3;   % width of stimulus bar (deg)
    settings.barheight       = .8;   % height of stimulus bar (deg)
    settings.ellipseArea     = settings.barwidth*settings.barheight; % ellipse size (deg^2)
    settings.jitter          = 0;   % amount of x/y-jitter (deg)
    settings.bgdac           = 128;  % background grayvalue (RGB)
    settings.fgdac           = 200;  % foreground grayvalue (RGB)
    settings.stimecc         = 7;    % stimulus eccentricity (deg)
    settings.ITT             = 1000;  % inter stimulus time (ms)

    % set bgdac
    settings.bgdac = 128;
    
    % Set feedback flag
    settings.feedback = 1;
    
    % screen info
    screenNumber=max(Screen('Screens'));       % use external screen if exists
    [w h]=Screen('WindowSize', screenNumber);  % screen resolution
    screen_resolution = [w h];                 % screen resolution
    screen_center = screen_resolution/2;       % screen center
    screen_distance = 60;                      % distance between observer and screen (in cm)
    screen_angle = 2*(180/pi)*(atan((settings.screen_width/2) / screen_distance)) ; % total visual angle of screen
    screen_ppd = screen_resolution(1) / screen_angle;  % pixels per degree
    screen_fixposxy = screen_resolution .* [.5 .5]; % fixation position
    settings.ellipseArea = .3;
    settings.ellipseArea = settings.ellipseArea * screen_ppd^2;
    
    % open screen
    gray=GrayIndex(screenNumber);
    windowPtr = Screen('OpenWindow',screenNumber,gray,[],32,2);
    
    % Set break points
    breakpoints = round((1:(breaknum-1)).* (settings.nTrials/breaknum));
    
    % Fixation information
    fixsize = 4;
    fixcol = 0;
    nextFlipTime = 0; % just to initialize...
    aborted = 0;
    
    % Create trials
    % 1-> target presence (1,0)
    % 2-> response (1,0)
    % 3-> real stim time (ms?)
    % 4-> reaction time (ms)
    % 5-> number of stimuli [setsize] (int)
    % 6-> delta
    % 7:14-> first array orientations (radians)
    % 15:22-> second array orientations (radians)
    % 23:30-> x locations (pixels)
    % 31:38-> y locations (pixels)
    % 39-> eccentricity of ellipses
       
    % If want to draw stimuli randomly using generative model
    % This is the default, how it's done in paper
    if(strcmp(drawstimuli,'randomly'))
        TrialMat = zeros(nTrials,39);
        % Set target presence/absence
        TrialMat(:,1) = rand(nTrials,1)>.5;
        % Set setsize
        TrialMat(:,5) = setsizeval(ceil(rand(nTrials,1)*length(setsizeval)));
        
        % Set deltas
        if (strcmp(deltatype,'constants'))
            TrialMat(:,6) = deltaval(ceil(rand(nTrials,1)*length(deltaval)));
        elseif (strcmp(deltavartype,'Uniform'))
            TrialMat(:,6) = (deltavartyperange*rand(nTrials,1));
        end
        
        % Set reliabilities
        TrialMat(:,39) = reliabilityval(ceil(rand(nTrials,1)*length(reliabilityval)));
        % Set location of change
        TrialMat(:,1) = ceil(rand(nTrials,1).*TrialMat(:,1).*TrialMat(:,5));
        
        % Generate stimuli orientations, uniform over theta (random)
        for i = 1:nTrials
            % Initial orientations
            randOrt = (pi)*rand(1,TrialMat(i,5));
            TrialMat(i,7:(TrialMat(i,5)+7-1)) = randOrt;
            
            % Second orientations
            % Initially the same as the first
            TrialMat(i,15:(TrialMat(i,5)+15-1)) = randOrt;
            
            % Add change at target location
            TrialMat(i,14+TrialMat(i,1))=TrialMat(i,14+TrialMat(i,1))+TrialMat(i,6)*(TrialMat(i,1)~=0);
            
            % Pick angle of first one rand and add the rest in counterclockwise
            % fashion with angle spacing = 2pi/max(setsizeval)
            locAngles = rand()*2*pi+(1:(max(setsizeval)))*(2*pi)/max(setsizeval);
            [X Y] = pol2cart(locAngles,screen_ppd * settings.stimecc);
            TrialMat(i,23:(22+TrialMat(i,5))) = X(1:TrialMat(i,5)) + screen_center(1);
            TrialMat(i,31:(30+TrialMat(i,5))) = Y(1:TrialMat(i,5)) + screen_center(2);
            
        end
        
    end
    
    % Jitter stimulus positions
    TrialMat(:,23:38) = TrialMat(:,23:38) + round((rand(nTrials,16)-.5)*settings.jitter*screen_ppd);
    
    % Create all stimulus patches, with degreeres degrees between each
    % possible stimulus orientations
    clear StimPatches;
    
    % Set the number of possible orientations based on resolution
    res = degreeres;
    numDegrees = ceil(180/res);
    
    % StimPatches holds the ellipse image patches
    StimPatches = zeros(length(reliabilityval),numDegrees);
    
    % create stimulus patches
    clear stimtex;           
    
    % Fill StimPatches
    radius = round(sqrt(settings.stimsize/pi)*screen_ppd);
    StimSizes = zeros(length(reliabilityval),numDegrees+1,2);
    colMat = getColMat();
    for i=1:length(reliabilityval)
        % Draw a patch for each orientation
        for j = 0:numDegrees

            template = fspecial('disk',radius);
            for jj=1:3
                template(template~=0)=colMat(round(j*res+1),jj);
                im(:,:,jj)=template;
            end
            im(im==0)=settings.bgdac;
            
            StimPatches(i,j+1) = Screen('MakeTexture',windowPtr,im);
            StimSizes(i,j+1,:) = size(im(:,:,1));
            StimSizes(i,j+1,:) = StimSizes(i,j+1,[2 1]);
            
        end
    end
        
    % Add on the actual degree values & StimPatch used for orientation
    % [exact_deg_i exact_deg_f StimPatch_idx_i StimPatch_dx_f]
    % Columns 56-63 and 64-71 contain the indices for getting images
    ExactDegs = [TrialMat(:,7:14)*(180/pi) TrialMat(:,15:22)*(180/pi)];
    ApproxDegs = [round(ExactDegs(:,1:8)*(1/res)) round(ExactDegs(:,9:16)*(1/res))];
    ApproxDegs = [mod(ApproxDegs(:,1:8),180)+1 mod(ApproxDegs(:,9:16),180)+1];
    TrialMat = [TrialMat ExactDegs ApproxDegs];
       
    % Randomize the trials
    TrialMat = TrialMat(randperm(nTrials),:);
    datafilename = ['output/' settings.subjid '_' settings.expID '_' 'color' '_' settings.disttype '_' num2str(settings.nTrials) '_' datestr(now,'yyyymmddTHHMMSS') '.mat'];
    
    % show start screen
    Screen('TextSize',windowPtr,20);
    textx = 100;
    texty = screen_center(2) - 50;
    dy = 30;
    Screen('DrawText',windowPtr,'Your task is to detect whether there is any change between',textx,texty,[255 255 255]); texty = texty + dy;
    Screen('DrawText',windowPtr,'the colors of the first and second sets of stimuli',textx,texty,[255 255 255]); texty = texty + 2*dy;
    Screen('DrawText',windowPtr,'The change can be of any magnitude. Press "y" if you think',textx,texty,[255 255 255]); texty = texty + dy;
    Screen('DrawText',windowPtr,'there is a change and "u" if you think they are the same.',textx,texty,[255 255 255]); texty = texty + 4*dy;   
    
    Screen('Flip', windowPtr);
    waitForKey;
    Screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    Screen('flip',windowPtr,nextFlipTime);
    nextFlipTime = GetSecs + 1;
    Screen('TextSize', windowPtr, 15);
    HideCursor;

    % Begin Trials!
    i = 0;
    while (i < nTrials) && ~aborted
        
        % If demo, then first few trials with long stimulus times
        if ~isempty(strfind(subjid,'DEMO')) && i<5
            stimtime = .5;
        else
            stimtime = .1;
        end        
        
        i = i+1;
        
        % SCREEN 1a: FIXATION
        Screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
        
        % First array, for each element
        for j= 1:TrialMat(i,5)
            squeeze(StimSizes(find(reliabilityval == TrialMat(i,39)),TrialMat(i,55+j),:));
            srcrect = [0 0 squeeze(StimSizes(find(reliabilityval == TrialMat(i,39)),TrialMat(i,55+j),:))'];
            destrect = CenterRectOnPoint(srcrect,TrialMat(i,22+j),TrialMat(i,30+j));
            Screen('drawtexture',windowPtr,StimPatches(find(reliabilityval == TrialMat(i,39)),TrialMat(i,55+j)),srcrect,destrect,0);
        end
        
        %Screen('DrawText',windowPtr,['Trial # = ' num2str(trialnr)],20,0,round(settings.bgdac*1.2));
        Screen('flip',windowPtr,nextFlipTime);
        nextFlipTime = GetSecs + stimtime;
        
        % SCREEN 1b: BLANK
        
        Screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
        Screen('flip',windowPtr,nextFlipTime);
        nextFlipTime = GetSecs + 1;
        
        % SCREEN 2: STIMULUS
        Screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
        
        % Second array, for each element
        for j= 1:TrialMat(i,5)
            srcrect = [0 0 squeeze(StimSizes(find(reliabilityval == TrialMat(i,39)),TrialMat(i,63+j),:))'];
            destrect = CenterRectOnPoint(srcrect,TrialMat(i,22+j),TrialMat(i,30+j));
            Screen('drawtexture',windowPtr,StimPatches(find(reliabilityval == TrialMat(i,39)),TrialMat(i,63+j)),srcrect,destrect,0);
        end
        
        Screen('flip',windowPtr,nextFlipTime);
        stimStartTime = GetSecs;
        nextFlipTime = GetSecs + stimtime; %settings.stimtime/1000;
        
        if (settings.makeScreenShot)
            grabSize = 2.5 * screen_ppd * settings.stimecc;
            grabrect = CenterRectOnPoint([0 0 grabSize grabSize],screen_fixposxy(1),screen_fixposxy(2));
            im = Screen('getimage',windowPtr,grabrect);
            imwrite(im,['screenshots/stim_' num2str(i) '.png'],'png');
        end
        
        % SCREEN 3a: YES/NO RESPONSE
        Screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
        Screen('flip',windowPtr,nextFlipTime);
        REALSTIMTIME = round(1000*(GetSecs - stimStartTime));
        
        yesKey = KbName('y');
        noKey = KbName('u');
        escKey = KbName('esc');
        responseStartTime = GetSecs;
        done=0;
        while ~done
            keyCode = waitForKey;
            if (keyCode==yesKey)
                YESNORESP = 1;
                done=1;
            elseif (keyCode==noKey)
                YESNORESP = 0;
                done=1;
            elseif (keyCode == escKey)
                aborted=1;
                break;
            end
        end
        
        CORRECT = (TrialMat(i,1) && keyCode==yesKey) || (~TrialMat(i,1) && keyCode==noKey);
        RT = round(1000*(GetSecs - responseStartTime));
        
        TrialMat(i,2) = YESNORESP;
        TrialMat(i,3) = REALSTIMTIME;
        TrialMat(i,4) = RT;
        
        save(datafilename,'settings','TrialMat');
        
        % SCREEN 4: INTER TRIAL DISPLAY
        Screen('fillRect',windowPtr,settings.bgdac);
        if (settings.feedback)
            if (CORRECT)
                drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),[0 200 0],fixsize);
            else
                drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),[200 0 0],fixsize);
            end
        else
            drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
        end
        Screen('flip',windowPtr);
        nextFlipTime = GetSecs + (settings.ITT/1000);
        
        % show progress info + target reminder
        
        if ~isempty(intersect(breakpoints,i))
            Screen('fillRect',windowPtr,settings.bgdac);
            Screen('flip',windowPtr,nextFlipTime);
            
            if (intersect(breakpoints,i) == breakpoints(ceil(breaknum/2)))
                TempTrialMat = TrialMat(1:i,:);
                HalfPC = sum((TempTrialMat(:,1)>0)-TempTrialMat(:,2)==0)/i;
                
                Screen('fillRect',windowPtr,settings.bgdac);
                Screen('DrawText',windowPtr,['You are halfway. You got ' num2str(HalfPC*100) '% correct so far. Press <ENTER> to continue'],250,screen_center(2) - 50,[255 255 255]);
                Screen('Flip', windowPtr);
                waitForKey;
                
                drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
                Screen('Flip', windowPtr);
                waitSecs(1.2);
                nextFlipTime = GetSecs + .5;
                
            else
                Screen('fillRect',windowPtr,settings.bgdac);
                tic
                
                cdtime = 60;   % time in seconds
                
                while toc<cdtime                    
                    Screen('DrawText',windowPtr,['You have finished ' num2str(round(100*i/settings.nTrials)) '% of the trials'],100,screen_center(2)-80,[255 255 255]);
                    Screen('DrawText',windowPtr,['Please take a short break now. You can continue in ' num2str(round(cdtime-toc)) ' seconds…'],100,screen_center(2)-100,[255 255 255]);                    
                    Screen('Flip', windowPtr);
                end                
                Screen('DrawText',windowPtr,'Press any key to continue',100,screen_center(2)-80,[255 255 255]);
                Screen('Flip', windowPtr);
                waitForKey;

                drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
                Screen('Flip', windowPtr);
                waitSecs(1.2);
                nextFlipTime = GetSecs + .5;
            end
            
        end
        
    end
    
    % Calculate percent correct
    PC = sum((TrialMat(1:i,1)>0)-TrialMat(1:i,2)==0)/i;
    
    % show start screen
    Screen('fillRect',windowPtr,settings.bgdac);
    Screen('DrawText',windowPtr,['End of this session. You got ' num2str(PC*100) '% correct. Press <ENTER> to continue'],250,screen_center(2) - 50,[255 255 255]);
    Screen('Flip', windowPtr);
    key = 0;
    while (key ~= 13)
        key = waitForKey;
    end
    
    %-%-%-%-%-%-%-
    %- FINALIZE %-
    %-%-%-%-%-%-%-
    ShowCursor;
    Screen('closeall');
    clear all_color;
    save all_color;
    
catch
    % catch error
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    clear all_color;
    save all_color;
    psychrethrow(psychlasterror);
    
end % try ... catch %


%-%-%-%-%-%-%-%-%-%-%-%-%- HELPER FUNCTIONS %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

function keyCode = waitForKey
keyCode = ones(1,256);
while sum(keyCode(1:254))>0
    [keyIsDown,secs,keyCode] = KbCheck;
end
while sum(keyCode(1:254))==0
    [keyIsDown,secs,keyCode] = KbCheck;
end
keyCode = min(find(keyCode==1));

function drawfixation(windowPtr,x,y,color,size)
Screen('DrawLine',windowPtr,color,x-size,y,x+size,y,2);
Screen('DrawLine',windowPtr,color,x,y-size,x,y+size,2);


function colMat = getColMat()

% This conversion is based on Bruce Lindbloom's excellent website about color spaces
% http://brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

% define stuff
ref_white  = [0.95047 1 1.08883];   % D65 reference white (this is a point in XYZ-space)
CIE_center = [58 12 13];   % (L,a,b) center point in CIE Lab space
CIE_radius = 50; % radius of the circle along which we'll pick our colors
epsilon = 216/24389; % some junction point in some function; see CIE documentation 
kappa = 24389/27; % another variable in the same function; see CIE documentation

% compute the (L,a,b)-values of the colorwheel 
L = CIE_center(1) * ones(1,181);
a = CIE_center(2) + CIE_radius*cos((0:180)*2*pi/180);
b = CIE_center(3) + CIE_radius*sin((0:180)*2*pi/180);

%-%-%-%-%-%-%-%-%-%-%-
% Step 1: Lab to XYZ %
%-%-%-%-%-%-%-%-%-%-%-
f_y = (L+16)/116;
f_x = a/500 + f_y;
f_z = f_y - b/200;

x_r                     = f_x.^3;
x_r(x_r<=epsilon)       = (116*f_x(x_r<=epsilon)-16)/kappa;
y_r(L>kappa*epsilon)    = ((L(L>kappa*epsilon)+16)/116).^3;
y_r(L<=kappa*epsilon)   = L(L<=kappa*epsilon)/kappa;
z_r                     = f_z.^3;
z_r(z_r<=epsilon)       = (116*f_z(z_r<=epsilon)-16)/kappa;

X = x_r * ref_white(1);
Y = y_r * ref_white(2);
Z = z_r * ref_white(3);

%-%-%-%-%-%-%-%-%-%-%-
% Step 2: XYZ to RGB %
%-%-%-%-%-%-%-%-%-%-%-

% From Bruce Lindbloom's website: "Given an XYZ color whose components are 
% in the nominal range [0.0, 1.0] and whose reference white is the same as 
% that of the RGB system, the% conversion to companded RGB is done in two steps."
%
% --> we need to check if our monitor's white point is D65 (that's the whitepoint that Wen-Chuang was using and that i'm using heret oo)

% STEP A: XYZ to Linear RGB
M = [ 3.2404542 -1.5371385 -0.4985314; ...   % conversion matrix for sRGB with D65 reference white 
     -0.9692660  1.8760108  0.0415560; ...   % source: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
      0.0556434 -0.2040259  1.0572252];
RGB = M*[X; Y; Z];
  
% STEP B: Companding
r = RGB(1,:); g = RGB(2,:); b = RGB(3,:);
R(r<=0.0031308) = 12.92*r(r<=0.0031308);
R(r >0.0031308) = 1.055*r(r >0.0031308).^(1/2.4);
G(g<=0.0031308) = 12.92*g(g<=0.0031308);
G(g >0.0031308) = 1.055*g(g >0.0031308).^(1/2.4);
B(b<=0.0031308) = 12.92*b(b<=0.0031308);
B(b >0.0031308) = 1.055*b(b >0.0031308).^(1/2.4);

colMat = round(255*[R; G; B]');
