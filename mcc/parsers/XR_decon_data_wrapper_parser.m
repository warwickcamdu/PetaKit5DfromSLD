function [] = XR_decon_data_wrapper_parser(dataPaths, varargin)


%#function XR_cudaDeconFrame3D
%#function XR_cppDeconFrame3D
%#function XR_RLdeconFrame3D
%#function XR_RotateFrame3D

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('deconPathstr', '',  @(x) ischar(x));
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 2) && islogical(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('Channels', [488, 560, 642], @(x) isnumeric(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', true, @(x) islogical(x) || ischar(x));
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false, @(x) numel(x) == 1 && islogical(x) || ischar(x));
ip.addParameter('onlyFirstTP', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Decon', true, @(x) islogical(x) || ischar(x));
ip.addParameter('cudaDecon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('cppDecon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('cppDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv', @ischar);
ip.addParameter('loadModules', 'module load gcc/4.8.5; module load fftw/3.3.6-gcc; module load boost/1.65.1-gcc; module load libtiff/4.1.0; ', @ischar);
ip.addParameter('cudaDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/cudaDeconv' , @ischar);
ip.addParameter('OTFGENPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/radialft' , @ischar); % point to radialft file
ip.addParameter('Background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dzPSF', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('EdgeErosion', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('ErodeByFTP', true, @(x) islogical(x) || ischar(x)); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('deconRotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('psfFullpaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('rotatePSF', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DeconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('OTFCumThresh', 0.9, @(x) isnumeric(x) || ischar(x)); % OTF cumutative sum threshold
ip.addParameter('hanWinBounds', [0.8, 1.0], @(x) isnumeric(x) || ischar(x)); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x)); % decon in skewed space
ip.addParameter('fixIter', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('errThresh', [], @(x) isnumeric(x) || ischar(x)); % error threshold for simplified code
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x)); % debug mode for simplified code
ip.addParameter('saveStep', 5, @(x) isnumeric(x) || ischar(x)); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @(x) islogical(x) || ischar(x)); % psf generation
ip.addParameter('GPUJob', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('BatchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % block size 
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeMethod', 'inmemory', @ischar); % inmemory, inplace. 
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('damper', 1, @(x) isnumeric(x) || ischar(x)); % damp factor for decon result
ip.addParameter('scaleFactor', [], @(x) isnumeric(x) || ischar(x)); % scale factor for decon result
ip.addParameter('deconOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset for decon result
ip.addParameter('deconMaskFns', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('maxWaitLoopNum', 10, @(x) isnumeric(x) || ischar(x)); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
deconPathstr = pr.deconPathstr;
Overwrite = pr.Overwrite;
ChannelPatterns = pr.ChannelPatterns;
Channels = pr.Channels;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
Save16bit = pr.Save16bit;
onlyFirstTP = pr.onlyFirstTP;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
Decon = pr.Decon;
cudaDecon = pr.cudaDecon;
cppDecon = pr.cppDecon;
cppDeconPath = pr.cppDeconPath;
loadModules = pr.loadModules;
cudaDeconPath = pr.cudaDeconPath;
OTFGENPath = pr.OTFGENPath;
Background = pr.Background;
dzPSF = pr.dzPSF;
EdgeErosion = pr.EdgeErosion;
ErodeByFTP = pr.ErodeByFTP;
deconRotate = pr.deconRotate;
psfFullpaths = pr.psfFullpaths;
rotatePSF = pr.rotatePSF;
DeconIter = pr.DeconIter;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
hanWinBounds = pr.hanWinBounds;
skewed = pr.skewed;
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;
GPUJob = pr.GPUJob;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
damper = pr.damper;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
deconMaskFns = pr.deconMaskFns;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
maxWaitLoopNum = pr.maxWaitLoopNum;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
GPUConfigFile = pr.GPUConfigFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(ChannelPatterns) && ~isempty(ChannelPatterns) && strcmp(ChannelPatterns(1), '{')
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(Channels)
    Channels = str2num(Channels);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(ObjectiveScan)
    ObjectiveScan = str2num(ObjectiveScan);
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = str2num(sCMOSCameraFlip);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(onlyFirstTP)
    onlyFirstTP = str2num(onlyFirstTP);
end
if ischar(parseSettingFile)
    parseSettingFile = str2num(parseSettingFile);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(Decon)
    Decon = str2num(Decon);
end
if ischar(cudaDecon)
    cudaDecon = str2num(cudaDecon);
end
if ischar(cppDecon)
    cppDecon = str2num(cppDecon);
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(dzPSF)
    dzPSF = str2num(dzPSF);
end
if ischar(EdgeErosion)
    EdgeErosion = str2num(EdgeErosion);
end
if ischar(ErodeByFTP)
    ErodeByFTP = str2num(ErodeByFTP);
end
if ischar(deconRotate)
    deconRotate = str2num(deconRotate);
end
if ischar(psfFullpaths) && ~isempty(psfFullpaths) && strcmp(psfFullpaths(1), '{')
    psfFullpaths = eval(psfFullpaths);
end
if ischar(rotatePSF)
    rotatePSF = str2num(rotatePSF);
end
if ischar(DeconIter)
    DeconIter = str2num(DeconIter);
end
if ischar(wienerAlpha)
    wienerAlpha = str2num(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(hanWinBounds)
    hanWinBounds = str2num(hanWinBounds);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(fixIter)
    fixIter = str2num(fixIter);
end
if ischar(errThresh)
    errThresh = str2num(errThresh);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(saveStep)
    saveStep = str2num(saveStep);
end
if ischar(psfGen)
    psfGen = str2num(psfGen);
end
if ischar(GPUJob)
    GPUJob = str2num(GPUJob);
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(damper)
    damper = str2num(damper);
end
if ischar(scaleFactor)
    scaleFactor = str2num(scaleFactor);
end
if ischar(deconOffset)
    deconOffset = str2num(deconOffset);
end
if ischar(deconMaskFns) && ~isempty(deconMaskFns) && strcmp(deconMaskFns(1), '{')
    deconMaskFns = eval(deconMaskFns);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(parseParfor)
    parseParfor = str2num(parseParfor);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(maxTrialNum)
    maxTrialNum = str2num(maxTrialNum);
end
if ischar(unitWaitTime)
    unitWaitTime = str2num(unitWaitTime);
end
if ischar(maxWaitLoopNum)
    maxWaitLoopNum = str2num(maxWaitLoopNum);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_decon_data_wrapper(dataPaths, deconPathstr=deconPathstr, Overwrite=Overwrite, ...
    ChannelPatterns=ChannelPatterns, Channels=Channels, SkewAngle=SkewAngle, ...
    dz=dz, xyPixelSize=xyPixelSize, Reverse=Reverse, ObjectiveScan=ObjectiveScan, ...
    sCMOSCameraFlip=sCMOSCameraFlip, Save16bit=Save16bit, onlyFirstTP=onlyFirstTP, ...
    parseSettingFile=parseSettingFile, flipZstack=flipZstack, Decon=Decon, ...
    cudaDecon=cudaDecon, cppDecon=cppDecon, cppDeconPath=cppDeconPath, loadModules=loadModules, ...
    cudaDeconPath=cudaDeconPath, OTFGENPath=OTFGENPath, Background=Background, ...
    dzPSF=dzPSF, EdgeErosion=EdgeErosion, ErodeByFTP=ErodeByFTP, deconRotate=deconRotate, ...
    psfFullpaths=psfFullpaths, rotatePSF=rotatePSF, DeconIter=DeconIter, RLMethod=RLMethod, ...
    wienerAlpha=wienerAlpha, OTFCumThresh=OTFCumThresh, hanWinBounds=hanWinBounds, ...
    skewed=skewed, fixIter=fixIter, errThresh=errThresh, debug=debug, saveStep=saveStep, ...
    psfGen=psfGen, GPUJob=GPUJob, BatchSize=BatchSize, BlockSize=BlockSize, ...
    zarrSubSize=zarrSubSize, largeFile=largeFile, largeMethod=largeMethod, ...
    zarrFile=zarrFile, saveZarr=saveZarr, damper=damper, scaleFactor=scaleFactor, ...
    deconOffset=deconOffset, deconMaskFns=deconMaskFns, parseCluster=parseCluster, ...
    parseParfor=parseParfor, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, maxWaitLoopNum=maxWaitLoopNum, ...
    mccMode=mccMode, ConfigFile=ConfigFile, GPUConfigFile=GPUConfigFile);

end
