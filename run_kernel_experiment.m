function output=run_kernel_experiment(mode,outputDir,checkpointDir)
if nargin<1,mode='full';end
if nargin<2,outputDir=fullfile(fileparts(mfilename('fullpath')),'results','section_6_3');end
if nargin<3,checkpointDir=fullfile(fileparts(mfilename('fullpath')),'checkpoints','section_6_3');end
if ~isfolder(outputDir),mkdir(outputDir);end
if ~isfolder(checkpointDir),mkdir(checkpointDir);end

definitions=case_definitions();
if strcmp(mode,'quick'),definitions=definitions(1);end
nRefine=max([definitions.level])-1;
seed=20260821;
hierarchy=fem_build_unstructured_red_hierarchy(4,10,nRefine,seed);
mg=mg_setup(hierarchy,[],struct('nSmooth',3,'tau',1.7,'useGalerkin',true));
[parts,kernelProperties]=cp_p16_representative_kernels();
summary=table;
pcgAll=table;
spectrumAll=table;

for ic=1:numel(definitions)
    def=definitions(ic);
    checkpoint=fullfile(checkpointDir,sprintf('case_%02d.mat',def.id));
    if isfile(checkpoint)
        saved=load(checkpoint,'caseSummary','casePcg','caseSpectrum');
        caseSummary=saved.caseSummary;
        casePcg=saved.casePcg;
        caseSpectrum=saved.caseSpectrum;
    else
        base=hierarchy.mesh{def.level};
        delta=def.ratio*base.hmax;
        fprintf('\nKernel case %d: N_h=%d, delta/h=%.4g\n', ...
            def.id,numel(base.free),def.ratio);
        ext=cp_extend_base_mesh(base,delta,max(4,ceil(1/base.hmax)));
        bundle=cp_assemble_kernel_bundle(ext,parts,delta);
        relativeVolumeError=abs(bundle{1}.stats.pairVolume- ...
            cp_p16_pair_volume(delta))/cp_p16_pair_volume(delta);
        assert(relativeVolumeError<2e-7);
        rows=struct([]);
        casePcg=table;
        caseSpectrum=table;
        for ik=1:numel(parts)
            nl=bundle{ik};
            nl.name=sprintf('case_%02d_%s',def.id,parts{ik}.kernel);
            opts=struct('level',def.level,'tol',1e-8,'maxit',300, ...
                'spectralLimit',1400,'directLimit',6000, ...
                'highModeLimit',2500,'timingRepeats',3, ...
                'practicalOnly',true,'seed',seed);
            result=mg_validate_paircell_case(nl,hierarchy,mg,opts);
            result.pcg.case_id(:)=def.id;
            result.spectrum.case_id(:)=def.id;
            casePcg=[casePcg;result.pcg]; %#ok<AGROW>
            caseSpectrum=[caseSpectrum;result.spectrum]; %#ok<AGROW>
            row=struct('case_id',def.id,'kernel',string(parts{ik}.kernel), ...
                'n_dof',numel(base.free),'h',base.hmax,'delta',delta, ...
                'delta_over_h',def.ratio, ...
                'lambda_min',result.spectrum.lambda_min, ...
                'lambda_max',result.spectrum.lambda_max, ...
                'condition_number',result.spectrum.condition_number, ...
                'maximum_PCG_iterations',max(result.pcg.pcg_iterations), ...
                'maximum_true_residual', ...
                max(result.pcg.recomputed_relative_residual));
            if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
        end
        caseSummary=struct2table(rows);
        save(checkpoint,'caseSummary','casePcg','caseSpectrum','def','-v7.3');
        clear bundle ext
    end
    summary=[summary;caseSummary]; %#ok<AGROW>
    pcgAll=[pcgAll;casePcg]; %#ok<AGROW>
    spectrumAll=[spectrumAll;caseSpectrum]; %#ok<AGROW>
end

paperTable=build_paper_table(summary);
writetable(paperTable,fullfile(outputDir,'table_6_4.csv'));
writetable(summary,fullfile(outputDir,'kernel_summary.csv'));
writetable(pcgAll,fullfile(outputDir,'pcg_details.csv'));
writetable(spectrumAll,fullfile(outputDir,'spectrum_details.csv'));
writetable(kernelProperties,fullfile(outputDir,'kernel_properties.csv'));
disp(paperTable);
output=struct('table',paperTable,'summary',summary, ...
    'pcg',pcgAll,'spectrum',spectrumAll,'kernelProperties',kernelProperties);
end

function definitions=case_definitions()
id=num2cell(1:7);
level=num2cell([2 3 3 3 4 4 4]);
ratio=num2cell([0.8 0.4 0.8 1.6 0.2 0.8 3.2]);
definitions=struct('id',id,'level',level,'ratio',ratio);
end

function paperTable=build_paper_table(summary)
names=["rho0","rhoS","rhoG","rho4","rhoC"];
prefix=["rho0","rho_step","rho_mix","rho4","rho_cone"];
caseIds=unique(summary.case_id,'stable');
paperTable=table;
for j=1:numel(caseIds)
    T=summary(summary.case_id==caseIds(j),:);
    row=table(T.n_dof(1),T.delta_over_h(1), ...
        'VariableNames',{'n_dof','delta_over_h'});
    for k=1:numel(names)
        S=T(T.kernel==names(k),:);
        row.([char(prefix(k)) '_condition_number'])=S.condition_number;
        row.([char(prefix(k)) '_maximum_PCG_iterations'])=S.maximum_PCG_iterations;
    end
    paperTable=[paperTable;row]; %#ok<AGROW>
end
end
