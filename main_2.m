function output=main_2(mode)
if nargin<1,mode='full';end
mode=validatestring(mode,{'full','quick'});
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
outputDir=fullfile(codeDir,'results','section_6_2');
checkpointDir=fullfile(codeDir,'checkpoints','section_6_2');
output=run_mesh_horizon_experiment(mode,outputDir,checkpointDir);
end
