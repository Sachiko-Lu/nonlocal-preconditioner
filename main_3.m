function output=main_3(mode)
if nargin<1,mode='full';end
mode=validatestring(mode,{'full','quick'});
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
outputDir=fullfile(codeDir,'results','section_6_3');
checkpointDir=fullfile(codeDir,'checkpoints','section_6_3');
output=run_kernel_experiment(mode,outputDir,checkpointDir);
end
