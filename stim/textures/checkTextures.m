% Check if modulating gabor has the same luminance at the different
% eccentricities with the different flip variants.

function checkTextures(texName,hemi,ecc)

% load textures
load(['./savedTextures/',texName,'.mat']);
myscreen.displayName = 'gdm';
initScreen(myscreen);

% loop through the different flip variants of the texture
for f = 1:length(textureParams.flip)
    mglClearScreen; mglFlush; mglClearScreen;
    tex = mglCreateTexture(texture{f,hemi,ecc,end});
    mglBltTexture(tex,[0 0]);
    mglFlush
    pause
end