% [texture, width, height] =
% make2ndOrderTex(carrierSF,gaborSF,gaborContrast,ecc,hemifield,target,varargin)
%
% varargin:
% screenName
% gaborAngle
% gaborPhase
% gaborSD
% width
% height
% dispFig

% written by Michael Jigo, October 2016
% This function creates 2nd-order textures using a first-order carrier and
% a Gabor modulator. This function will make a texture that is usable by
% MGL.

% The input arguments "ecc" and "hemifield" relate to the eccentricity and
% hemifield in which the gabor modulator will be displayed
% ecc is in visual degrees
% hemifield can be an integer between 1 and 2. hemifield=1 is the left
% hemifield; hemifield=2 is the right hemifield

% The "target" input argument determines whether a gabor will be displayed
% or not. I did this to try and fix luminance issue (i.e., modulated noise
% is usually of higher luminance than the carrier noise). It did not fix
% the issue.

function [texture, width, height] = make2ndOrderTex(carrierSF,gaborSF,gaborContrast,ecc,hemifield,target,varargin)
getArgs(varargin,{'screenName=gdm','gaborAngle=90','gaborPhase=0',...
    'gaborSD=0.8','width=30.5','height=10','dispFig=0','makeSmall=3',...
    'flipCarrier=0','newCarrier=1'});

%% Initialize parameters of stimulus
% initialize MGL to get screen parameters
if dispFig
    mglClose;
    initScreen(screenName);
end

% make fullscreen stimulus by default
if isempty(width) || isempty(height)
    width = mglGetParam('deviceWidth')-makeSmall;
    height = mglGetParam('deviceHeight')-makeSmall;
end

% calculate Cartesian coordinates of Gabor modulator position
switch hemifield
    case 2 % rightward along the horizontal meridian
        hemifield = 0;
    case 1 % leftward along the horizontal meridian
        hemifield = pi;
end
%% Make stimulus
% load saved carrier noise (this is done to ensure that the same noise is
% used in all trials)
if newCarrier==1
    carrier = makeCarrier(width,height,carrierSF);
else
    load(newCarrier); % loading a carrier based on its file name in newCarrier
end

% to avoid local adaptation, flip the carrier noise according to
% flipCarrier
switch flipCarrier
    case 1
        carrier = fliplr(carrier);
    case 2
        carrier = flipud(carrier);
    case 3
        carrier = fliplr(flipud(carrier));
end

% make preliminary gabor
% 1) Gaussian is made such that it extends to, at most, gaborSize
% 2) Equation taken from Methods of Sutter et al., 1995
grating = mglMakeGrating(width,height,gaborSF,gaborAngle,gaborPhase);
% gaussian window
window = mglMakeGaussian(width,height,gaborSD,gaborSD); % 1
% cosine window (doesn't really work)
% window = makeCosWindow(width,height,gaborSize,gaborSize);

% make gabor and a gabors with the highest and lowest possible contrast
gabor = gaborContrast*window.*grating;
gabor = 1+gabor; % 2

if ecc>0 || hemifield==0
    % move gabor
    gabor = moveGabor(gabor,ecc,hemifield);
end

% normalize the carrier relative to -2 and 2, the minimum and maximum
% of the stimuli
modCar = {carrier carrier.*gabor};
modCar = normZeroToOne(modCar,[-2 2]);
carrier = modCar{1}; modCar = modCar{2};
if target==1
    % make 2nd-order texture
    texture = 255*modCar;
elseif target==2
    % only make gabor
    texture = 255*normZeroToOne(gabor);
else
    texture = 255*carrier;
end

%% Display stimulus
if dispFig
    % display texture with mgl
    x = mglCreateTexture(texture);
    mglBltTexture(x,[0 0]);
    
    % for testing purposes, make rings at each eccentricity
    %     allEcc = [0 1.2 3.6 6 8.4];
    %             allEcc = [0 2.4 4.8 7.2 9.6];
    % uncomment if you want the center of the target position to be displayed
    % via a red dot
    %     [gaborXCent, gaborYCent] = pol2cart(hemifield,allEcc); % output is in dva
    %     mglFillOval(gaborXCent,gaborYCent,[0.25 0.25],[1 0 0]);
    
    % plot cues
    %     [x0, y0] = pol2cart(hemifield,0.35);
    %     [x1, y1] = pol2cart(hemifield,1);
    %     mglLines2(x0,y0,x1,y1,2.5,[0 0 0]);
    %     mglTextDraw(num2str(3),[0,0]);
    
    mglFlush;
    
    % disply carrier and modulated carrier, ensure that they're both the
    % same luminance
    %     figure('Name','modulated carrier');
    %     imagesc(texture); colormap gray; truesize; caxis([0 255]);
    %     title(sprintf('avg. luminance val. = %.2f',mean(texture(:))));
    %
    %     figure('Name','carrier');
    %     imagesc(255*carrier); colormap gray; truesize; caxis([0 255]);
    %     title(sprintf('avg. luminance val. = %.2f',mean(255*carrier(:))));
end

function movedGab = moveGabor(gab,ecc,theta)
% when mglMakeGaussian moves the center of the gaussian via xCenter and
% yCenter, the phase of the gabor patch changes. this results in different
% looking gabors presented when I move it to different locations.

% This code will simply translate the gabor to the desired position defined
% by eccentricty and theta.

% get cartesian coordinates of new location in visual angle
[x, y] = pol2cart(theta,ecc);

% convert to pixels
xDeg2Pix = mglGetParam('xDeviceToPixels');
yDeg2Pix = mglGetParam('yDeviceToPixels');
x = round(x*xDeg2Pix);
y = round(y*yDeg2Pix);
% equate pixel calculation with carrier
x = x+mod(x+1,2);
y = y+mod(y+1,2);

%% move Gaussian by adding the mode value
% I chose to use the mode because the majority of the image is not the
% gabor but gray space
nX = abs(x);
nY = abs(y);
addX = repmat(mode(gab(:)),size(gab,1),nX);
addY = repmat(mode(gab(:)),nY,size(gab,2));

% determine whether to add pixels to the left of right of image (X dim)
if x>0
    % move gaussian to right
    movedGab = [addX gab];
    movedGab(:,end:-1:end-nX+1) = [];
elseif x<0
    % move gaussian to left
    movedGab = [gab addX];
    movedGab(:,1:nX) = [];
end

% determine whether to add pixels to the top or bottom of image (Y dim)
if y>0
    % move gaussian up
    movedGab = [addY; movedGab];
    movedGab(end:-1:end-nY+1,:) = [];
elseif y<0
    % move gaussian down
    movedGab = [movedGab; addY];
    movedGab(1:nY,:) = [];
end
