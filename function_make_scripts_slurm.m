function function_make_scripts_slurm(cfg)

if ~isfield(cfg,'subi')
    cfg.subi=1;
end
if ~isfield(cfg,'subs')
    cfg.subs=[];
end
if ~isfield(cfg,'infile')
    cfg.infile=[];
end
if ~isfield(cfg,'nii')
    cfg.nii=[];
end
if ~isfield(cfg,'mask')
    cfg.mask=[];
end
if ~isfield(cfg,'regressor')
    cfg.regressor=[];
end
if ~isfield(cfg,'clusterstat')
    cfg.clusterstat='maxsum';
end
if ~isfield(cfg,'alpha')
    cfg.alpha=0.05;
end
if ~isfield(cfg,'cdtP')
    cfg.cdtP=0.01;
end
if ~isfield(cfg,'toi')
    cfg.toi=[];
end
if ~isfield(cfg,'outdir')
    cfg.outdir=[];
end


jobind=cfg.ind;
disp('Making jobs ...')

fid = fopen(['./jobs/job_' num2str(jobind) '.sh'],'w');

fprintf(fid,'#!/bin/bash\n\n');

fprintf(fid,['#SBATCH -p ' cfg.partition '\n']);
fprintf(fid,['#SBATCH -t ' cfg.time '\n']);
fprintf(fid,['#SBATCH -o ' './logs/log_' num2str(jobind) '\n']);
fprintf(fid,'#SBATCH --qos=normal\n');
% fprintf(fid,'#SBATCH --exclusive\n');
fprintf(fid,['#SBATCH --mem-per-cpu=' cfg.mem '\n\n']);

fprintf(fid,['matlab -nojvm -r "cd ' pwd '/; cfg.infile=''' cfg.infile '''; cfg.nii=''' cfg.nii '''; cfg.outdir=''' cfg.outdir '''; cfg.mask=''' cfg.mask '''; cfg.regressor=''' cfg.regressor '''; cfg.clusterstat=''' cfg.clusterstat '''; cfg.alpha=' num2str(cfg.alpha) '; cfg.cdtP=' num2str(cfg.cdtP) '; cfg.toi=[' num2str(cfg.toi) ']; cfg.subs=''' cfg.subs '''; cfg.subi=' num2str(cfg.subi) '; ' cfg.func ' ; exit;"']);

fclose(fid);

end
