% Written by Barbara Montagna, NYU, October 2007
% Edited by Michael Jigo to make it work with MGL, October 2016

% This program was written to create isotropic visual noise that is
% bandpass filtered in the spatial frequency domain. Visual noise is
% created by randomly sampling values from a normal distribution with mu=0
% and sigma=1. These random numbers are then associated with poins in the
% complex Fourier domain. A bandpass filter is applied, such that there is
% perfect transmission within a band of frequencies and zero transmission
% outside that band. The resulting noise contains an octave-wide band of
% spatial frequencies around a given central frequency. In the filtered
% noise, pixel values above or below 3 standard deviations from the mean
% are truncated and set to + or - 3 standard deviations from the mean,
% respectively. To achieve the largest possible contrast in teh stimulus,
% the noise images are normalized to produce values ranging from 0 to 1
% that will span the maximum range of the monitor used.

% Input:
% width and height of stimulus in visual degrees
% centralSF is the spatial frequency of the stimulus in cycles per degree
% (cpd)

% Output:
% outputs carrier in the same size as the gabor created by using
% mglMakeGrating and mglMakeGaussian

function carrier = makeCarrier(width,height,centralSF,dispFig)

%% Set parameters
if ieNotDefined('dispFig')
    dispFig = 0;
end
% in cycles per degree this is the only value to be changed to create noise 
% with different central frequencies
centralFrequency = centralSF; 

% visual degree to pixels
%xDeg2pix = mglGetParam('xDeviceToPixels');
%yDeg2pix = mglGetParam('yDeviceToPixels');

% TEMPORARY FOR MAKING FIGURES FOR POSTER
xDeg2pix = 1;
yDeg2pix = 1;

% pixels to visual degree
%xDegPerPix = mglGetParam('xPixelsToDevice');
%yDegPerPix = mglGetParam('yPixelsToDevice');

xDegPerPix = 1;
yDegPerPix = 1;
degPerPix = [xDegPerPix yDegPerPix];

% determine size of stimulus in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

% get a grid of x and y coordinates that has
% the correct number of pixels
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
[xMesh,~] = meshgrid(x,y);

%% FREQUENCY NOISE
% the initial noise image in the spatial domain
pixelsPerImage = size(xMesh); % in pixels (better if odd to avoid center to be shifted in the fourier transform)
noise = (rand(pixelsPerImage)*2)-1; % in the contrast domain, range of values: from -1 to 1

% central frequency and one octave bandwidth
cpdHigh = centralFrequency*(4/3);  % in cycles per degree
cpdLow = centralFrequency*(2/3); % in cycles per degree

% Conversions from cycles per degree to cycles per image
cyclesPerPixHigh = degPerPix*cpdHigh;
cyclesPerPixLow = degPerPix*cpdLow;
cyclesPerImageHigh = cyclesPerPixHigh.*pixelsPerImage;
cyclesPerImageLow = cyclesPerPixLow.*pixelsPerImage;

% Setting fLow=0 and fHigh=Inf will produce an all-pass filter.
fLow = cyclesPerImageLow;
fHigh =  cyclesPerImageHigh;
fNyquist = pixelsPerImage./2; % in cycles per image
filter = Bandpass2(pixelsPerImage,min(fLow./fNyquist),max(fHigh./fNyquist));
ft=filter.*(fftshift(fft2(noise)));
myNoise = ifft2(ifftshift(ft));

% truncate pixels that are above or below 3 SD from the mean
truncateCutoff = 3*std(myNoise(:));%3*std(reshape(myNoise,pixelsPerImage(1)*pixelsPerImage(2),1));
indexesPositiveVal = myNoise > truncateCutoff;
myNoise(indexesPositiveVal) = truncateCutoff;
indexesNegativeVal = myNoise < (-truncateCutoff);
myNoise(indexesNegativeVal) = -truncateCutoff;

% Normalize the Noise to achieve a highest possible contrast level
carrier =  myNoise./(max(abs(myNoise(:)))); % range from -1 to 1

% plot the initial noise in the spatial domain, the bandpass filter, the
% filtered noise in the fourier domain, and the filtered noise in the
% spatial domain
if dispFig
    figure(1);
    subplot(2,2,1); imagesc(noise); title('noise'); axis square
    subplot(2,2,2); imagesc(filter); title('filter'); colormap gray; axis square
    subplot(2,2,3); imagesc(abs(ft)); title('filtered noise in fourier domain'); axis square
    subplot(2,2,4); imagesc(carrier); title('filtered noise in the spatial domain'); axis square
    
    % plot the filtered noise on a scale from 0 to 255
    figure(2);
    subplot(1,2,1);surf(128+(carrier*128)); axis square
    subplot(1,2,2);imagesc(128+(carrier*128)); colormap gray; axis square
    
    % plot the filtered noise in its true size
    figure(3);
    imagesc(128+(carrier*128)); colormap gray; truesize
end
