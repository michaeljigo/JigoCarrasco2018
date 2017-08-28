% INPUTS:
% preMakeTextures(width,height,carrierSF,gaborSF,gaborPhase,gaborSD,hemifield,eccentricity,contrast)

% pre-make the textures to be used in the texture segmentation experiment

% hemifield, ecc, and contrast can be vectors. the function will make a
% texture at each value

function [texture, textureParams] = preMakeTextures(width,height,carrierSF,gaborSF,gaborPhase,gaborSD,gaborAngle,hemifield,eccentricity,contrast,flip,saveName)

% initialize screen
myscreen.displayName = 'gdm';
initScreen(myscreen);

% get the various carriers that will be used
carrierDir = ['./savedTextures/carrier',num2str(carrierSF),'CPD/'];
carriers = dir([carrierDir,'*.mat']);

% initialize variable
% flip = [0=original; 1=fliplr; 2=flipud; 3=flipud+fliplr]

% start counter
n = 0; disppercent(-inf);
nTotal = length(carriers)*length(flip)*length(hemifield)*length(eccentricity)*length(contrast);

% loop through each hemifield, ecc, and contrast value and make a texture
% for each combination one. Do this for every flipped variant of the
% carrier noise
for car = 1:length(carriers)
    texture = cell(length(flip),length(hemifield)+1,length(eccentricity),length(contrast));
    for f = 1:length(flip)
        for h = 1:length(hemifield)
            hemi = hemifield(h);
            for e = 1:length(eccentricity)
                ecc = eccentricity(e);
                for c = 1:length(contrast)
                    con = contrast(c);
                    temp = make2ndOrderTex(carrierSF,...
                        gaborSF,con,ecc,hemi,1,'width',width,'height',height,...
                        'gaborPhase',gaborPhase,'gaborSD',gaborSD,'flipCarrier',...
                        flip(f),'newCarrier',[carrierDir,carriers(car).name],...
                        'gaborAngle',gaborAngle);
                    
                    % make texture into format expected by openGL to hopefully
                    % speed things up (as instructed on MGL How Tos website)
                    temp = repmat(uint8(temp),[1 1 3]);
                    temp = permute(temp,[3 2 1]);
                    temp(4,:,:) = 255;
                    texture{f,h,e,c} = temp;
                    
                    % update counter
                    n = n+1;
                    disppercent(n/nTotal,sprintf('%i/%i textures made',n,nTotal));
                end
            end
        end
        
        % add target absent texture
        temp = make2ndOrderTex(carrierSF,gaborSF,...
            con,ecc,hemi,0,'width',width,'height',height,'gaborPhase',gaborPhase,...
            'gaborSD',gaborSD,'flipCarrier',flip(f),'newCarrier',[carrierDir,carriers(car).name]);
        temp = repmat(uint8(temp),[1 1 3]);
        temp = permute(temp,[3 2 1]);
        temp(4,:,:) = 255;
        texture{f,h+1,1,1} = temp; % add in noise as third hemifield
    end
    
    % store ecc, contrast, and hemifield so that I can make sure that I'm
    % loading the correct textures when loading into the experiment AND to
    % ensure that I'm indexing the correct textures
    textureParams.size = [width height];
    textureParams.carrierSF = carrierSF;
    textureParams.gaborSF = gaborSF;
    textureParams.gaborPhase = gaborPhase;
    textureParams.gaborSD = gaborSD;
    textureParams.hemifield = hemifield;
    textureParams.eccentricity = eccentricity;
    textureParams.contrast = contrast;
    textureParams.flip = flip;
    textureParams.gaborAngle = gaborAngle;
    textureParams.cellOrganization = {'flip' 'hemifield' 'eccentricity' 'contrast'};
    
    % save textureParams
    saveDir = './savedTextures/';
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    save([saveDir,saveName,num2str(car),'.mat'],'texture','textureParams','-v7.3');
end
disppercent(inf);