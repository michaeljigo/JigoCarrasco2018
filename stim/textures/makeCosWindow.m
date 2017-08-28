% Make a 2D cosine window. Used with mglMakeGrating to create a "top-hat"
% gabor.

% Input:
% width = width of image
% height = height of image
% winWidth = width of window
% winHeight = height of window

% All inputs are in visual degrees

% Modeled after mglMakeGaussian

function m = makeCosWindow(width,height,winWidth,winHeight)

% check arguments
m = [];
if nargin<4
  help makeCosWindow
  return
end

if ieNotDefined('xCenter'),xCenter = 0;end
if ieNotDefined('yCenter'),yCenter = 0;end

% defaults for xDeg2pix
if ieNotDefined('xDeg2pix')
  if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
  end
  xDeg2pix = mglGetParam('xDeviceToPixels');
end

% defaults for yDeg2pix
if ieNotDefined('yDeg2pix')
  if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
  end
  yDeg2pix = mglGetParam('yDeviceToPixels');
end

% get size of image in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

% get size of window in pixels
winWidthPixels = round(winWidth*xDeg2pix);
winHeightPixels = round(winHeight*yDeg2pix);
winWidthPixels = winWidthPixels + mod(winWidthPixels+1,2);
winHeightPixels = winHeightPixels + mod(winHeightPixels+1,2);

% get a grid of x and y coordinates that has 
% the correct number of pixels
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
[xMesh,yMesh] = meshgrid(x,y);

% now create a grid from the 1D cosine windows (using hann because it goes
% to zero)
x = blackman(winWidthPixels);
y = blackman(winHeightPixels);
[winX,winY] = meshgrid(x,y);

% center the window on the image
zeroX = find(xMesh(1,:)==0);
zeroY = find(yMesh(:,1)==0);
winIdxY = linspace(zeroY-floor(winHeightPixels/2),zeroY+floor(winHeightPixels/2),winHeightPixels);
winIdxX = linspace(zeroX-floor(winWidthPixels/2),zeroX+floor(winWidthPixels/2),winWidthPixels);
% verify that window will be centered on 0
if any(~isequal(median(winIdxY),zeroY)||~isequal(median(winIdxX),zeroX))
    error('window is not centered at (0,0)');
end

% prepare the image for the window
xMesh(winIdxY, winIdxX) = nan;
yMesh(winIdxY, winIdxX) = nan;

% make sure everything outside the window is 0
xMesh(~isnan(xMesh)) = 0;
yMesh(~isnan(yMesh)) = 0;

% put in the window
xMesh(winIdxY, winIdxX) = winX;
yMesh(winIdxY, winIdxX) = winY;

% finalize window
m = xMesh+yMesh;
