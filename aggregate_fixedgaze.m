addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130')

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

gazedistthresh = 20;

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

%% print genotype
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
genotypes = cell(size(nwbsessions));
tic
for ises = 1:numel(nwbsessions)
nwbfiles = cat(1, dir([datadir nwbsessions{ises} '/*.nwb']), dir([datadir nwbsessions{ises} '/*/*.nwb']));

% take filename  with shortest length or filename that does not contain probe
[~, fileind] = min(cellfun(@length, {nwbfiles.name}));
nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);

nwb = nwbRead(nwbspikefile);

genotypes{ises} = nwb.general_subject.genotype;
disp(nwb.general_subject.genotype)
end
toc

% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=2;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

allfields = fieldnames(ICsig_fixedgaze.ICwcfg1_presentations);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),2);
end
ICsigfields = allfields(validfields);
ICsigfieldsize2 = size2fields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};

allfields = fieldnames(RFCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCI_fixedgaze.(allfields{f}),2);
end
RFCIfields = allfields(validfields);
RFCIfieldsize2 = size2fields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};

allfields = fieldnames(RFCIspin_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCIspin_fixedgaze.(allfields{f}),2);
end
RFCIspinfields = allfields(validfields);
RFCIspinfieldsize2 = size2fields(validfields);

allfields = fieldnames(sizeCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(sizeCI_fixedgaze.(allfields{f}),2);
end
sizeCIfields = allfields(validfields);
sizeCIfieldsize2 = size2fields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};

allfields = fieldnames(oriparams_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(oriparams_fixedgaze.(allfields{f}),2);
end
oriparamsfields = allfields(validfields);
oriparamsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};

allfields = fieldnames(ori4params_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ori4params_fixedgaze.(allfields{f}),2);
end
ori4paramsfields = allfields(validfields);
ori4paramsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};


unit_amplitude_agg = cell(1,Nsessions);
unit_isi_violations_agg = cell(1,Nsessions);
unit_wfdur_agg = cell(1,Nsessions);
unit_amplitude_cutoff_agg = cell(1,Nsessions);
unit_presence_ratio_agg = cell(1,Nsessions);

neuoindagg = cell(numel(probes),Nsessions);
probeneuagg = cell(numel(probes),Nsessions);
neulocagg = cell(numel(probes),Nsessions);
neupeakchagg = cell(numel(probes),Nsessions);
neuctxagg = cell(numel(probes),Nsessions);
sesneuagg = cell(numel(probes),Nsessions);

meanFRvecagg = cell(size(probes));
sponFRvecagg = cell(size(probes));

for b = 1:numel(ICblocks)
    ICsig_fixedgazeagg.(ICblocks{b}) = cell(numel(probes), Nsessions);
end
RFCI_fixedgazeagg = cell(numel(probes), Nsessions);
RFCIspin_fixedgazeagg = cell(numel(probes), Nsessions);
sizeCI_fixedgazeagg = cell(numel(probes), Nsessions);
oriparams_fixedgazeagg = cell(numel(probes), Nsessions);
ori4params_fixedgazeagg = cell(numel(probes), Nsessions);

validfixedgazesessions = true(1,Nsessions);
for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    if ~exist([pathpp 'trackmouseeye.mat'], 'file')
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        validfixedgazesessions(ises) = false;
        continue
    end
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    load([pathpp 'qc_units.mat']) %'unit_amplitude', 'unit_isi_violations', 'unit_wfdur', 'unit_amplitude_cutoff', 'unit_presence_ratio'

    unit_amplitude_agg{ises} = unit_amplitude;
    unit_isi_violations_agg{ises} = unit_isi_violations;
    unit_wfdur_agg{ises} = unit_wfdur;
    unit_amplitude_cutoff_agg{ises} = unit_amplitude_cutoff;
    unit_presence_ratio_agg{ises} = unit_presence_ratio;
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        tic
%         load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'


            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        
        
%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        if nnz(floor(unit_peakch/1000)==probeind-1)==0
            fprintf('Probe %s Area %s: NO UNITS!!!\n', probes{iprobe}, visareas{iprobe} )
            continue
        end
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
%         disp(unique(probelocs)')
        
        load(sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        neuoindagg{iprobe, ises} = neuoind;
        probeneuagg{iprobe, ises} = iprobe*ones(length(neuoind),1);
        neulocagg{iprobe, ises} = neuloc;
        neupeakchagg{iprobe, ises} = unit_peakch(neuoind);
        neuctxagg{iprobe, ises} = neuctx;

        sesneuagg{iprobe, ises} = ises*ones(length(neuoind),1);
        % sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));

        meanFRvecagg{iprobe, ises} = meanFRvec';
        sponFRvecagg{iprobe, ises} = sponFRvec';
        
        Nneuoi = length(neuoind);
        
        for b = 1:numel(ICblocks)
            for f= 1:numel(ICsigfields)
                if ~isfield(ICsig_fixedgaze.(ICblocks{b}), ICsigfields{f})
                    tempmat = NaN( Nneuoi, ICsigfieldsize2(f) );
                else
                    tempmat = ICsig_fixedgaze.(ICblocks{b}).(ICsigfields{f});
                end
            ICsig_fixedgazeagg.(ICblocks{b}){iprobe,ises}.(ICsigfields{f}) = tempmat;
            end
        end
        
        % RFCIfields = fieldnames(RFCI);
        for f= 1:numel(RFCIfields)
            if ~isfield(RFCI_fixedgaze, RFCIfields{f})
                tempmat = NaN( Nneuoi, RFCIfieldsize2(f) );
            else
                tempmat = RFCI_fixedgaze.(RFCIfields{f});
            end
            RFCI_fixedgazeagg{iprobe,ises}.(RFCIfields{f}) = tempmat;
        end
        
        for f= 1:numel(RFCIspinfields)
            if ~isfield(RFCIspin_fixedgaze, RFCIspinfields{f})
                tempmat = NaN( Nneuoi, RFCIspinfieldsize2(f) );
            else
                tempmat = RFCIspin_fixedgaze.(RFCIspinfields{f});
            end
            RFCIspin_fixedgazeagg{iprobe,ises}.(RFCIspinfields{f}) = tempmat;
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        for f= 1:numel(sizeCIfields)
            if ~isfield(sizeCI_fixedgaze, sizeCIfields{f})
                tempmat = NaN( Nneuoi, sizeCIfieldsize2(f) );
            else
                tempmat = sizeCI_fixedgaze.(sizeCIfields{f});
            end
            sizeCI_fixedgazeagg{iprobe,ises}.(sizeCIfields{f}) = tempmat;
        end
        
        for f= 1:numel(oriparamsfields)
            if ~isfield(oriparams_fixedgaze, oriparamsfields{f})
                tempmat = NaN( Nneuoi, oriparamsfieldsize2(f) );
            else
                tempmat = oriparams_fixedgaze.(oriparamsfields{f});
            end
            oriparams_fixedgazeagg{iprobe,ises}.(oriparamsfields{f}) = tempmat;
        end
        
        for f= 1:numel(ori4paramsfields)
            if ~isfield(ori4params_fixedgaze, ori4paramsfields{f})
                tempmat = NaN( Nneuoi, ori4paramsfieldsize2(f) );
            else
                tempmat = ori4params_fixedgaze.(ori4paramsfields{f});
            end
            ori4params_fixedgazeagg{iprobe,ises}.(ori4paramsfields{f}) = tempmat;
        end
        
        toc
    end
    
end

Nrfs = RFCIfieldsize2(strcmp(RFCIfields, 'Rrfclassic'));
Nszs = sizeCIfieldsize2(strcmp(sizeCIfields, 'Rsizeclassic'));
Ndirs = oriparamsfieldsize2(strcmp(oriparamsfields, 'Rori'));
Noris = Ndirs/2;

ises=4;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

%%
unit_amplitude_all = cat(1,unit_amplitude_agg{:});
unit_isi_violations_all = cat(1,unit_isi_violations_agg{:});
unit_wfdur_all = cat(1,unit_wfdur_agg{:});
unit_amplitude_cutoff_all = cat(1,unit_amplitude_cutoff_agg{:});
unit_presence_ratio_all = cat(1,unit_presence_ratio_agg{:});

neuoindall = cat(1,neuoindagg{:});
probeneuall = cat(1,probeneuagg{:});
neulocall = cat(1,neulocagg{:});
neupeakchall = cat(1,neupeakchagg{:});
neuctxall = cat(1,neuctxagg{:});
sesneuall = cat(1,sesneuagg{:});
meanFRvecall = cat(1,meanFRvecagg{:});
sponFRvecall = cat(1,sponFRvecagg{:});
% figure; plot(sesneuall*10+probeneuall, '.')

ICsig_fixedgazeall = struct();
for b = 1:numel(ICblocks)
    temp = cat(1,ICsig_fixedgazeagg.(ICblocks{b}){:});
    for f= 1:numel(ICsigfields)
        ICsig_fixedgazeall.(ICblocks{b}).(ICsigfields{f}) = cat(1,temp.(ICsigfields{f}));
    end
end

RFCI_fixedgazeall = struct();
temp = cat(1,RFCI_fixedgazeagg{:});
for f= 1:numel(RFCIfields)
    RFCI_fixedgazeall.(RFCIfields{f}) = cat(1,temp.(RFCIfields{f}));
end

RFCIspin_fixedgazeall = struct();
temp = cat(1,RFCIspin_fixedgazeagg{:});
for f= 1:numel(RFCIspinfields)
    RFCIspin_fixedgazeall.(RFCIspinfields{f}) = cat(1,temp.(RFCIspinfields{f}));
end

sizeCI_fixedgazeall = struct();
temp = cat(1,sizeCI_fixedgazeagg{:});
for f= 1:numel(sizeCIfields)
    sizeCI_fixedgazeall.(sizeCIfields{f}) = cat(1,temp.(sizeCIfields{f}));
end

oriparams_fixedgazeall = struct();
temp = cat(1,oriparams_fixedgazeagg{:});
for f= 1:numel(oriparamsfields)
    oriparams_fixedgazeall.(oriparamsfields{f}) = cat(1,temp.(oriparamsfields{f}));
end

ori4params_fixedgazeall = struct();
temp = cat(1,ori4params_fixedgazeagg{:});
for f= 1:numel(ori4paramsfields)
    ori4params_fixedgazeall.(ori4paramsfields{f}) = cat(1,temp.(ori4paramsfields{f}));
end

clearvars temp

%% sanity check: neuoind is ordered according to the probe order
for ises = 1:Nsessions
    if validfixedgazesessions(ises)
    if ~isequal( cat(1,neuoindagg{:,ises}), (1:length(cat(1,neuoindagg{:,ises})))' )
        error('neuoind is not ordered according to probe indices')
    end
    end
end

if ~isequal(floor(neupeakchall/1000), probeneuall-1)
    error('check probeneuall')
end

%%
save([datadir 'postprocessed\openscope_popavg_fixedgazeall_', num2str(gazedistthresh), 'pix.mat'], ...
    'gazedistthresh', 'validfixedgazesessions', 'probes', 'visblocks', 'ICblocks', ...
    'nwbsessions', 'neuoindall', 'probeneuall', 'neulocall', 'neupeakchall', 'sesneuall', 'neuctxall', ...
    'meanFRvecall', 'sponFRvecall', 'vis', ...
    'unit_amplitude_all', 'unit_isi_violations_all', 'unit_wfdur_all', 'unit_amplitude_cutoff_all', 'unit_presence_ratio_all', ...
    'ICsigfields', 'ICsigfieldsize2', 'ICsig_fixedgazeall', ...
    'RFCIfields', 'RFCIfieldsize2', 'RFCIfieldsize2', 'RFCI_fixedgazeall', ...
    'RFCIspinfields', 'RFCIspinfieldsize2', 'RFCIspin_fixedgazeall', ...
    'sizeCIfields', 'sizeCIfieldsize2', 'sizeCI_fixedgazeall', ...
    'oriparamsfields', 'oriparamsfieldsize2', 'oriparams_fixedgazeall', ...
    'ori4paramsfields', 'ori4paramsfieldsize2', 'ori4params_fixedgazeall', '-v7.3')

%% aggregate Ron Roff and psthagg
aggpsth = true;
if ~aggpsth
probesR = {'C'};
else
    probesR = probes;
end

load([pathpp 'postprocessed_probeC.mat'], 'psthtli', 'vis')
Nneuronsall = numel(sesneuall);

vistrialtypes_fixedgazeagg = struct();
vistrialrep_fixedgazeagg = struct();
vistrialorder_fixedgazeagg = struct();
vistrial_fixedgazeagg = struct();
if aggpsth
psthavg_fixedgazeall = struct();
% psthsem_fixedgazeall = struct();
end
Ronavg_fixedgazeall = struct();
Roffavg_fixedgazeall = struct();
Ronstd_fixedgazeall = struct();
Roffstd_fixedgazeall = struct();
for b = 1:numel(visblocks)
    temptt = unique(vis.(visblocks{b}).trialorder);
    if strcmp(visblocks{b}, 'RFCI_presentations')
        temptt = unique(vis.(visblocks{b}).trialorder(1:4:end));
    end
    Ntt = numel(temptt);
    for iprobe = 1:numel(probesR)
        if aggpsth
            psthavg_fixedgazeall.(visblocks{b}) = NaN(length(psthtli), Ntt, Nneuronsall);
            % psthsem_fixedgazeall.(visblocks{b}) = NaN(length(psthtli), Ntt, Nneuronsall);
        end
        Ronavg_fixedgazeall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Roffavg_fixedgazeall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Ronstd_fixedgazeall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Roffstd_fixedgazeall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
    end
end

skippedsessions = false(Nsessions,1);
for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    try
        load([pathpp 'trackmouseeye.mat'])
    catch
        skippedsessions(ises) = true;
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end

    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probesR)
        tic
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        if ~ismember(probesR{iprobe}, probelist.probes)
            warning('%s: Probe%s does not exist', nwbsessions{ises}, probesR{iprobe})
            continue
        end

        tempneusesprobe = sesneuall==ises & probeneuall==probeind;

        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probesR{iprobe}))
        % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        %     load(sprintf('%svisresponses_probe%s.mat', pathpp, probesR{iprobe}))
        %     % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
                
        if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
            error('check neuoind')
        end
        Nneuoi = length(neuoind);

        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        
        fprintf('Probe %s Area %s: %d/%d\n', probesR{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
        disp(unique(probelocs)')
        
%         probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
%         neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
%         neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
%         neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
%         
%         sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
%         sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
        
        
        
        % visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        %     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
        %ICblocks: each stim is 0.4s, inter-trial interval is 0.4s, static images
        %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
        for b = 1:numel(visblocks)
            vistrialtypes_fixedgazeagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
            vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            if strcmp(visblocks{b}, 'RFCI_presentations')
                vistrialtypes_fixedgazeagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder(1:4:end));
                vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder(1:4:end);
            end
            vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            
            if ismember(visblocks{b}, ICblocks)
                tlon = psthtli>0 & psthtli<=400;
                tloff = psthtli>400 & psthtli<=800;
                temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
            elseif strcmp(visblocks{b}, 'RFCI_presentations')
                tlon = psthtli>0 & psthtli<=1000;
                tloff = psthtli>1000 & psthtli<=1000;
                tempgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
                templikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
                temptrialsfixedgaze = tempgazedist<gazedistthresh & ~templikelyblink;
                temptrialsfixedgaze = reshape(repmat(temptrialsfixedgaze', 4,1),[],1);
            elseif strcmp(visblocks{b}, 'sizeCI_presentations')
                tlon = psthtli>0 & psthtli<=250;
                tloff = psthtli>250 & psthtli<=750;
                temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
            else
                error('vis block not recognized')
            end
            
            vistrial_fixedgazeagg(ises).(visblocks{b}) = temptrialsfixedgaze;
            Ntt = numel(vistrialtypes_fixedgazeagg(ises).(visblocks{b}));
            vistrialrep_fixedgazeagg(ises).(visblocks{b}) = zeros(Ntt,1);

            temppsth = NaN(length(psthtli), Ntt, Nneuoi);
            tempRonavg = NaN(Nneuoi, Ntt);
            tempRoffavg = NaN(Nneuoi, Ntt);

            temppsthsem = NaN(length(psthtli), Ntt, Nneuoi);
            tempRonstd = NaN(Nneuoi, Ntt);
            tempRoffstd = NaN(Nneuoi, Ntt);
            for ii = 1:Ntt
                trialsoi = temptrialsfixedgaze & vis.(visblocks{b}).trialorder==vistrialtypes_fixedgazeagg(ises).(visblocks{b})(ii);
                vistrialrep_fixedgazeagg(ises).(visblocks{b})(ii) = nnz(trialsoi);

                temppsth(:,ii,:) = mean(1000*psth.(visblocks{b})(:,trialsoi,:), 2);
                tempRonavg(:,ii) = squeeze(mean(1000*psth.(visblocks{b})(tlon,trialsoi,:), [1 2]));
                tempRoffavg(:,ii) = squeeze(mean(1000*psth.(visblocks{b})(tloff,trialsoi,:), [1 2]));

                temppsthsem(:,ii,:) = std(1000*psth.(visblocks{b})(:,trialsoi,:), 0,2)/sqrt(nnz(trialsoi));
                tempRonstd(:,ii) = squeeze(std( mean(1000*psth.(visblocks{b})(tlon,trialsoi,:),1), 0,2));
                tempRoffstd(:,ii) = squeeze(std( mean(1000*psth.(visblocks{b})(tloff,trialsoi,:),1), 0,2));
            end

            if aggpsth
                psthavg_fixedgazeall.(visblocks{b})(:,:,tempneusesprobe) = temppsth;
                % psthsem_fixedgazeall.(visblocks{b})(:,:,tempneusesprobe) = temppsthsem;
            end
            Ronavg_fixedgazeall.(visblocks{b})(tempneusesprobe,:) = tempRonavg ;
            Roffavg_fixedgazeall.(visblocks{b})(tempneusesprobe,:) = tempRoffavg ;
            Ronstd_fixedgazeall.(visblocks{b})(tempneusesprobe,:) = tempRonstd ;
            Roffstd_fixedgazeall.(visblocks{b})(tempneusesprobe,:) = tempRoffstd ;
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        toc
    end
    
end

for b = 1:numel(visblocks)
    tic
    invalidsessions1 = unique(sesneuall( squeeze( all(isnan(psthavg_fixedgazeall.(visblocks{b})), [1,2]) ) ));
    % invalidsessions2 = unique(sesneuall( squeeze( all(isnan(psthsem_fixedgazeall.(visblocks{b})), [1,2]) ) ));
    invalidsessions2 = invalidsessions1;
    invalidsessions3 = unique(sesneuall( all(isnan(Ronavg_fixedgazeall.(visblocks{b})), 2) ));
    invalidsessions4 = unique(sesneuall( all(isnan(Roffavg_fixedgazeall.(visblocks{b})), 2) ));
    invalidsessions5 = unique(sesneuall( all(isnan(Ronstd_fixedgazeall.(visblocks{b})), 2) ));
    invalidsessions6 = unique(sesneuall( all(isnan(Roffstd_fixedgazeall.(visblocks{b})), 2) ));

    if strcmp(visblocks{b}, 'RFCI_presentations') && isequal(invalidsessions1, invalidsessions2, invalidsessions3, invalidsessions5)
        disp([visblocks{b} ' following sessions are invalid'])
        disp(invalidsessions1)
    elseif isequal(invalidsessions1, invalidsessions2, invalidsessions3, invalidsessions4, invalidsessions5, invalidsessions6)
        disp([visblocks{b} ' following sessions are invalid'])
        disp(invalidsessions1)
    else
        error('%s: check invalid sessions', visblocks{b})
    end
    toc
end

if aggpsth
    save([datadir 'postprocessed\openscope_psthavg_fixedgazeall_', num2str(gazedistthresh), 'pix.mat'], ...
        'validfixedgazesessions', 'probes', 'visareas', 'visind', 'nwbsessions', ...
        'neuoindall', 'probeneuall', 'neulocall', 'neupeakchall', 'sesneuall', 'neuctxall', ...
        'unit_amplitude_all', 'unit_isi_violations_all', 'unit_wfdur_all', 'unit_amplitude_cutoff_all', 'unit_presence_ratio_all', ...
        'vistrialtypes_fixedgazeagg', 'vistrialrep_fixedgazeagg', 'vistrialorder_fixedgazeagg', 'vistrial_fixedgazeagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavg_fixedgazeall', 'Ronavg_fixedgazeall', 'Roffavg_fixedgazeall', ...
        'Ronstd_fixedgazeall', 'Roffstd_fixedgazeall', '-v7.3') %'psthsem_fixedgazeall', 
end

%% MATCH BETWEEN CRF POSITION FOR ALL TRIALS VS FIXED GAZE TRIALS: ~63% neurons show match
popavg = load([datadir 'postprocessed\openscope_popavg_all.mat']);

neuinarea = strcmp(neulocall, 'VISp2/3'); % 86% neurons show match
neuvalfg = ismember(popavg.sesneuall, sesneuall);
RFindclassic = popavg.RFCIall.RFindclassic(neuvalfg);

%neuoi = neuinarea;
neuoi = neuinarea & popavg.RFCIall.Pkw_rfclassic(neuvalfg)<0.05;
figure; histogram2(RFindclassic(neuoi), RFCI_fixedgazeall.RFindclassic(neuoi), 'displaystyle', 'tile')

disp('match between CRF position for all trials vs fixed gaze trials')
disp(mean(RFindclassic(neuoi)==RFCI_fixedgazeall.RFindclassic(neuoi)))

%%
copyfile('S:\OpenScopeData\00248_v240130\postprocessed\openscope_popavg_fixedgazeall_20pix.mat', 'G:\My Drive\RESEARCH\ICexpts_revision23\')
copyfile('S:\OpenScopeData\00248_v240130\postprocessed\openscope_psthavg_fixedgazeall_20pix.mat', 'G:\My Drive\RESEARCH\ICexpts_revision23\')
