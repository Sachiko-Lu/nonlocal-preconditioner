function output=run_mesh_horizon_experiment(mode,outputDir,checkpointDir)
if nargin<1,mode='full';end
if nargin<2,outputDir=fullfile(fileparts(mfilename('fullpath')),'results','section_6_2');end
if nargin<3,checkpointDir=fullfile(fileparts(mfilename('fullpath')),'checkpoints','section_6_2');end
if ~isfolder(outputDir),mkdir(outputDir);end
if ~isfolder(checkpointDir),mkdir(checkpointDir);end

if strcmp(mode,'quick')
    nRefine=2;
    cases=case_definitions(false);
else
    nRefine=5;
    cases=case_definitions(true);
end

seed=20260821;
hierarchy=fem_build_unstructured_red_hierarchy(4,10,nRefine,seed);
mg=mg_setup(hierarchy,[],struct('nSmooth',3,'tau',1.7,'useGalerkin',true));
pcgAll=table;
spectrumAll=table;
uniqueRows=struct([]);

for j=1:numel(cases)
    def=cases(j);
    checkpoint=fullfile(checkpointDir,[def.name,'.mat']);
    if isfile(checkpoint)
        saved=load(checkpoint,'result','row');
        result=saved.result;
        row=saved.row;
    else
        base=hierarchy.mesh{def.level};
        delta=def.ratio*base.hmax;
        fprintf('\n%s: N_h=%d, delta/h=%.4g\n',def.name,numel(base.free),def.ratio);
        [A,stats]=cp_assemble_p16_constant_fast(base,delta,2e-11);
        nl=struct('Adbl',A,'delta',delta,'beta',1/8,'name',def.name, ...
            'kernel','rho0','m',16,'q',0.5,'geometryCertified',false);
        if def.level<=5
            opts=struct('level',def.level,'tol',1e-8,'maxit',300, ...
                'spectralLimit',5000,'directLimit',6000, ...
                'highModeLimit',5000,'timingRepeats',3, ...
                'practicalOnly',~def.compare,'seed',seed);
            result=mg_validate_paircell_case(nl,hierarchy,mg,opts);
        else
            result=validate_large_practical_case(A,hierarchy,mg,delta,def.level,def.name);
        end
        practical=result.spectrum(string(result.spectrum.preconditioner)== ...
            "symmetric_MG_plus_true_lump",:);
        pcg=result.pcg(string(result.pcg.preconditioner)== ...
            "symmetric_MG_plus_true_lump",:);
        row=struct('case_name',string(def.name),'path',string(def.path), ...
            'level',def.level,'n_dof',numel(base.free),'h',base.hmax, ...
            'delta',delta,'delta_over_h',def.ratio, ...
            'lambda_min',practical.lambda_min,'lambda_max',practical.lambda_max, ...
            'condition_number',practical.condition_number, ...
            'maximum_PCG_iterations',max(pcg.pcg_iterations), ...
            'maximum_true_residual',max(pcg.recomputed_relative_residual), ...
            'relative_pair_volume_error',stats.relativePairVolumeError);
        save(checkpoint,'result','row','stats','def','-v7.3');
        clear A
    end
    result.pcg.case_name(:)=string(def.name);
    result.spectrum.case_name(:)=string(def.name);
    pcgAll=[pcgAll;result.pcg]; %#ok<AGROW>
    spectrumAll=[spectrumAll;result.spectrum]; %#ok<AGROW>
    if isempty(uniqueRows),uniqueRows=row;else,uniqueRows(end+1)=row;end %#ok<AGROW>
end

uniqueTable=struct2table(uniqueRows);
summary=expand_shared_case(uniqueTable);
comparison=preconditioner_comparison(pcgAll,spectrumAll);
writetable(summary,fullfile(outputDir,'table_6_2.csv'));
writetable(comparison,fullfile(outputDir,'table_6_3.csv'));
writetable(pcgAll,fullfile(outputDir,'pcg_details.csv'));
writetable(spectrumAll,fullfile(outputDir,'spectrum_details.csv'));
plot_mesh_horizon_results(summary,outputDir);
disp(summary(:,{'path','n_dof','h','delta_over_h','condition_number', ...
    'maximum_PCG_iterations','maximum_true_residual'}));
disp(comparison);
output=struct('summary',summary,'comparison',comparison, ...
    'pcg',pcgAll,'spectrum',spectrumAll);
end

function cases=case_definitions(fullRun)
cases=struct('name','shared','path','shared','level',2,'ratio',0.8,'compare',true);
if fullRun,lastLevel=6;else,lastLevel=3;end
for ell=3:lastLevel
    cases(end+1)=struct('name',sprintf('constant_L%d',ell-1), ...
        'path','constant','level',ell,'ratio',0.8,'compare',false); %#ok<AGROW>
    cases(end+1)=struct('name',sprintf('increasing_L%d',ell-1), ...
        'path','increasing','level',ell,'ratio',0.8*2^(ell-2),'compare',false); %#ok<AGROW>
    cases(end+1)=struct('name',sprintf('decreasing_L%d',ell-1), ...
        'path','decreasing','level',ell,'ratio',0.8/2^(ell-2),'compare',false); %#ok<AGROW>
end
end

function summary=expand_shared_case(uniqueTable)
shared=uniqueTable(uniqueTable.path=="shared",:);
rows=table;
paths=["constant","increasing","decreasing"];
for k=1:numel(paths)
    first=shared;
    first.path(:)=paths(k);
    rest=uniqueTable(uniqueTable.path==paths(k),:);
    rows=[rows;first;rest]; %#ok<AGROW>
end
summary=rows;
end

function comparison=preconditioner_comparison(pcg,spectrum)
keep=string(spectrum.case_name)=="shared";
spectrum=spectrum(keep,:);
names=["K_inverse","ideal_K_plus_consistent_mass", ...
    "symmetric_MG_plus_true_lump"];
labels=["K_h^{-1}","K_h^{-1}+beta*delta^2*M_h^{-1}","B_h"];
rows=struct([]);
for j=1:numel(names)
    s=spectrum(string(spectrum.preconditioner)==names(j),:);
    p=pcg(string(pcg.case_name)=="shared" & ...
        string(pcg.preconditioner)==names(j),:);
    row=struct('preconditioner',labels(j), ...
        'condition_number',s.condition_number, ...
        'minimum_PCG_iterations',min(p.pcg_iterations), ...
        'maximum_PCG_iterations',max(p.pcg_iterations));
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
comparison=struct2table(rows);
end
