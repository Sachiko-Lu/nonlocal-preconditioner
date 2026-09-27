function results=main_1(mode)
if nargin<1,mode="full";end
mode=validatestring(mode,{'full','quick'});
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
outputDir=fullfile(codeDir,'results','section_6_1');
if strcmp(mode,'quick')
    nValues=[4 8];
else
    nValues=[4 8 16 32 64];
end
results=run_convergence_experiment(nValues,outputDir,false);
end
