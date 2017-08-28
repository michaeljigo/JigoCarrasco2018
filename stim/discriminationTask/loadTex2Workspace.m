function workspaceTex = loadTex2Workspace(expName)
% clear workspace
clearvars -except workspaceTex expName

% check if the textures are loaded
if ieNotDefined('workspaceTex')
    % load the texture into the workspace
    texDir = ['./textures/savedTextures/',expName,'/'];
    tex = dir([texDir,'*.mat']);
    if isempty(tex)
        error([expName,' is not a valid experiment name.']);
    end
    % because I'm using different carriers and due to saving/loading times of
    % large files, I am going to load each "piece" of the total texture
    % structure and then concatenate them
    for t = 1:length(tex)
        fprintf('Loading texture: %i/%i \n',t,length(tex));
        load([texDir,tex(t).name]);
        % concatenate
        startIdx = (t-1)*length(textureParams.flip)+1;
        endIdx = startIdx+length(textureParams.flip)-1;
        allTex(startIdx:endIdx,:,:,:) = texture;
    end
    workspaceTex.textures = allTex;
    workspaceTex.texParams = textureParams;
    
    clearvars -except workspaceTex
    fprintf('Done.\n');
end