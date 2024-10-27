% PERHAPS CCG flank should have been calculated with jitter-corrected CCG?

%{
node strength can be calculated with:
-CRF edge potential
-Pearson correlation
-Okun's population coupling
-mutual information between pairs of neurons

sensory sensitivity can be cauculated with:
-CRF AUC
-SP_ICvsRC
-mutual information between vis stim and neural responses
%}
% IMPORTANT NOT TO DO GENPATH
% addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
ses2anal = [1 6];
for ises = ses2anal
    clearvars -except datadir nwbsessions ses2anal ises
    %disp(nwbsessions{ises})
    fprintf('%d %s\n', ises, nwbsessions{ises})
    sesclk = tic;
    
    neuopt = 'ctx';
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed_probeC.mat'], 'vis')
    
    pathsv = [datadir 'CCG' filesep nwbsessions{ises} filesep];
    if ~exist(pathsv, 'dir')
        mkdir(pathsv)
    end
    
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location'
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data'
%     load([pathpp 'spiketimes.mat'])
    load([pathpp 'visresponses.mat'])
    % isequal(ststartend, [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1])
    Tres = 0.001; % 1ms
    
    neuctx = contains(neuallloc, 'VIS');
    neuctxind = find(neuctx);
    
    Nneuctx = numel(neuctxind);
    neulocctx = neuallloc(neuctxind);
    % isequal(neulocctx, find(contains(neuloc, 'VIS')))
    
    %% extract spike times
    if exist([pathpp 'spiketimes.mat'], 'file')
        load([pathpp 'spiketimes.mat'])
    else
    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    nwb = nwbRead(nwbspikefile); %, 'ignorecache');
    unit_times_data = nwb.units.spike_times.data.load();
    
    recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
    visarealabels = zeros(Nneuctx,1);
    for a = 1:numel(recarealabels)
        if strcmp(recarealabels{a}, 'VISp')
            neuinarea = contains(neulocctx, 'VISp') & ~contains(neulocctx, 'VISpm');
        else
            neuinarea = contains(neulocctx, recarealabels{a});
        end
        visarealabels(neuinarea) = a;
    end
    
    ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
    ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
    ICblockstartend = floor(ICblockstartend/Tres)+1;
    
    Nneurons = length(unit_ids);
    spiketimes = cell(Nneurons, 1);
    last_idx = 0;
    for ii = 1:Nneurons
        unit_id = unit_ids(ii);
        
        %     assert(unit_trials_idx(i) == unit_times_idx(i), 'Expected unit boundaries to match between trials & spike_times jagged arrays')
        start_idx = last_idx + 1;
        end_idx = unit_times_idx(ii);
        
        spiketimes{ii} = unit_times_data(start_idx:end_idx);
        
        last_idx = end_idx;
    end
    
    Tres = 0.001; % 1ms
    stlen = ceil((max(unit_times_data)+1)/Tres); % add 1s buffer/padding after the last spike timing
    
    disp([stlen, Nneurons])
    save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuallloc', 'neuctx', 'neuctxind', 'neulocctx', ...
        'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')
    end

    %%    
    switch neuopt
        case 'ctx'
            neuctx = contains(neuallloc, 'VIS');
        case 'ctrctx'
            neuctx = contains(neuallloc, 'VIS') & RFCIall.RFindclassic==1 & RFCIall.Pkw_rfclassic<0.05;
        case 'ctxL23'
            neuctx = contains(neuallloc, 'VIS') & contains(neuallloc, '2/3');
    end
    neuctxind = find(neuctx) ;
    
    %%
    ctxspiketrain = false(Nneuctx, ststartend(end) );
    for ii = 1:Nneuctx
        ci = neuctxind(ii);
        ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
    end
    
    ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));
    % ctxspiketrain = ctxspiketrain(:,ststartend(1):ststartend(2));
    
    stlen = size(ctxspiketrain, 2);
    spkcntvec = sum(ctxspiketrain,2);
    sqrtspkcntmat = sqrt( spkcntvec * spkcntvec' );

    
    %%
    NFFT = 2^ceil(log2(stlen));
    CCGffttl = -NFFT/2:NFFT/2-1;
    Thalfwin = 100;
    CCGtli_fft = -Thalfwin:Thalfwin;
    fprintf('Number of neurons: %d\n', Nneuctx)
    
    tic
    ctxCCG_fft = NaN(Nneuctx, Nneuctx, length(CCGtli_fft));
    for ci = 1:Nneuctx
        for i100 = 1:ceil((Nneuctx-ci+1)/100)
            tempneuvec = ci+100*(i100-1): min([Nneuctx ci-1+100*(i100)]);
            %disp([tempneuvec(1) tempneuvec(end)])
            tempCCGfft = fftshift(ifft(fft(ctxspiketrain(tempneuvec,:)', NFFT) .* conj(fft(ctxspiketrain(ci,:)', NFFT))),1);
            ctxCCG_fft(ci,tempneuvec,:) = tempCCGfft(ismember(CCGffttl, CCGtli_fft),:)';
            ctxCCG_fft(tempneuvec,ci,:) = tempCCGfft(ismember(CCGffttl, CCGtli_fft),:)';
        end
        if mod(ci,10)==1
            toc
            fprintf('done with neuron #%d\n', ci)
        end
    end
    ctxCCGfft = ctxCCG_fft./sqrtspkcntmat;
    toc
    disp('calculated ctxCCG')
    
    switch neuopt
        case 'ctx'
            save([pathsv 'ctxCCGfft.mat'], 'stlen', 'spkcntvec', 'neuctx', 'neulocctx', 'visarealabels', 'CCGtli_fft', 'ctxCCGfft', '-v7.3')
        case 'ctrctx'
            neuctrctx = neuctx;
            neulocctrctx = neulocctx;
            ctrctxCCGfft = ctxCCGfft;
            ctrctxCCGweight = ctxCCGweight;
            save([pathsv 'ctrctxCCGfft.mat'], 'stlen', 'spkcntvec', 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli_fft', 'ctrctxCCGfft', '-v7.3')
        case 'ctxL23'
            neuctxL23 = neuctx;
            neulocctxL23 = neulocctx;
            ctxL23CCGfft = ctxCCGfft;
            ctxL23CCGweight = ctxCCGweight;
            save([pathsv 'ctxL23CCGfft.mat'], 'stlen', 'spkcntvec', 'neuctxL23', 'neulocctxL23', 'visarealabels', 'CCGtli_fft', 'ctxL23CCGfft', '-v7.3')
    end
    
    
    toc(sesclk)
end

%% match sessions between this dataset and v240130 and move corresponding sessions
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsesv240130 = {nwbdir.name};
nwbsesv240130 = nwbsesv240130(~contains(nwbsesv240130, 'Placeholder') & ...
    ( contains(nwbsesv240130, 'sub-') | contains(nwbsesv240130, 'sub_') ));
neulocaggv240130 = cell(size(nwbsesv240130));
for ises = 1:numel(nwbsesv240130)
    clearvars neuloc unit_peakch electrode_location electrode_id
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed\' nwbsesv240130{ises} '\'];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    neuloc = electrode_location(revmapelecid(unit_peakch+1));
    neulocaggv240130{ises} = neuloc;
end

datadir = 'S:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsesv0 = {nwbdir.name};
nwbsesv0 = nwbsesv0(~contains(nwbsesv0, 'Placeholder') & ...
    ( contains(nwbsesv0, 'sub-') | contains(nwbsesv0, 'sub_') ));
neulocaggv0 = cell(size(nwbsesv0));
for ises = 1:numel(nwbsesv0)
    clearvars neuloc unit_peakch electrode_location electrode_id
    pathpp = ['S:\OpenScopeData\000248\postprocessed\' nwbsesv0{ises} '\'];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    electrode_location = cellstr(electrode_location);
    neuloc = electrode_location(revmapelecid(unit_peakch+1));
    neulocaggv0{ises} = neuloc;
end

sesindsv240130 = zeros(size(nwbsesv240130));
for ises = 1:numel(nwbsesv240130)
    c = cellfun(@isequal, neulocaggv0, repmat(neulocaggv240130(ises),size(neulocaggv0)), 'UniformOutput',false);
    sesv0ind = find(cat(1,c{:}));
    if numel(sesv0ind)==1
        sesindsv240130(ises) = sesv0ind;
    else
        warning('%d %s does not have a direct match with the old version', ises, nwbsesv240130{ises})
    end
end
% Warning: 1 sub-619293 does not have a direct match with the old version 
% Warning: 6 sub-625554 does not have a direct match with the old version 

for ises = 1:numel(nwbsesv240130)
    if sesindsv240130(ises)==0
        warning('%d %s does not have a direct match with the old version', ises, nwbsesv240130{ises})
        continue
    end
    pathccg = ['S:\OpenScopeData\000248\CCG\' nwbsesv0{sesindsv240130(ises)} filesep];
    ccgfnv0 = [pathccg 'ctxCCGfft.mat'];
    % save([pathsv 'ctxCCGfft.mat'], 'stlen', 'spkcntvec', 'neuctx', 'neulocctx', 'visarealabels', 'CCGtli_fft', 'ctxCCGfft', '-v7.3')
    pathsv = ['S:\OpenScopeData\00248_v240130\CCG' filesep nwbsesv240130{ises} filesep];
    mkdir(pathsv)
    ccgfnv240130 = [pathsv 'ctxCCGfft.mat'];
    copyfile(ccgfnv0, pathsv)
end

whos('-file', 'S:\OpenScopeData\00248_v240130\CCG\sub-619293\ctxCCGfft.mat')
whos('-file', 'S:\OpenScopeData\00248_v240130\CCG\sub-637484\ctxCCGfft.mat')

%% approximate ctxCCGsm0 from ctxCCGfft by averaging -12~12ms window
load('S:\OpenScopeData\000248\CCG\sub_1175512783\ctxCCGsm0.mat')
load('S:\OpenScopeData\000248\CCG\sub_1175512783\ctxCCGfft.mat')

% note, high correlation
figure; hold all
plot(ctxCCGsm0, mean(ctxCCGfft(:,:,CCGtli_fft>=-12 & CCGtli_fft<=12),3), 'o')
xl=xlim;
plot(xl,xl,'k-')

%% NOTE ERROR IN CALCULATING ctxCCGfft -- when row>column, CCG vector should be flipped!!!
load('S:\OpenScopeData\000248\CCG\sub_1175512783\ctxCCGfft.mat')
load('S:\OpenScopeData\000248\CCG\sub_1175512783\ctxCCG.mat')

isequaln(ctxCCG, ctxCCGfft(:,:,ismember(CCGtli_fft, CCGtli))) % not true

r=107; c= 121; % test for arbitrary pair
isequal(ctxCCG(r,c,:), flip(ctxCCG(c,r,:))) % true

N = size(ctxCCGfft,1);
ctxCCGfftnew = ctxCCGfft;
tic
for t = 1:length(CCGtli_fft)
    tempCCG = ctxCCGfft(:,:,t);
    tempCCGflip = ctxCCGfft(:,:,length(CCGtli_fft)+1-t);
    tempCCG(tril(true(N),-1)) = tempCCGflip(tril(true(N),-1));
    ctxCCGfftnew(:,:,t) = tempCCG;
end
toc

isequaln(ctxCCG, ctxCCGfftnew(:,:,ismember(CCGtli_fft, CCGtli))) % true if not for floating point errors
max(abs(ctxCCG - ctxCCGfftnew(:,:,ismember(CCGtli_fft, CCGtli))),[],'all')
figure; hold all
plot(reshape(ctxCCG,[],1), reshape(ctxCCGfftnew(:,:,ismember(CCGtli_fft, CCGtli)),[],1), '.')
xl=xlim; plot(xl,xl,'k-')

r=107; c= 121; % test for arbitrary pair
isequal(ctxCCGfftnew(r,c,:), flip(ctxCCGfftnew(c,r,:))) % true

figure; hold all
plot(CCGtli_fft, squeeze(ctxCCGfftnew(r,c,:)), 'k-')
plot(CCGtli, squeeze(ctxCCG(r,c,:)), 'r--')


datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions
    tic
    pathccg = ['S:\OpenScopeData\00248_v240130\CCG\' nwbsessions{ises} filesep];
    load([pathccg 'ctxCCGfft.mat'])
    N = size(ctxCCGfft,1);
    ctxCCGfftnew = ctxCCGfft;
    for t = 1:length(CCGtli_fft)
        tempCCG = ctxCCGfft(:,:,t);
        tempCCGflip = ctxCCGfft(:,:,length(CCGtli_fft)+1-t);
        tempCCG(tril(true(N),-1)) = tempCCGflip(tril(true(N),-1));
        ctxCCGfftnew(:,:,t) = tempCCG;
    end
    save([pathccg 'ctxCCGfftnew.mat'], 'stlen', 'spkcntvec', 'neuctx', 'neulocctx', 'visarealabels', 'CCGtli_fft', 'ctxCCGfftnew', '-v7.3')
    toc
end
