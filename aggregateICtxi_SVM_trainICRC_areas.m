addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

neuopt = 'RS';
svmdesc = 'trainICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

pathsv = [datadir 'SVM_' svmdesc '_selectareas' filesep];
loadallsvm = true;
allblocks = true;

visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
if allblocks
    ICblocknames = {'ICkcfg0','ICkcfg1','ICwcfg0','ICwcfg1'};
else
    ICblocknames = {'ICwcfg1'};
end

% % initialize
whichvisarea = visareas{1};
whichICblock = ICblocknames{1};
pathsvm = [pathsv nwbsessions{1} filesep];
svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
load(svmfn, 'SVMtrainICRC')
traintrialtypes = reshape( SVMtrainICRC.trialtypes, 1,[]);
Ntt = numel(traintrialtypes);
Nsplits=size(SVMtrainICRC.spkcnt.traintrialinds,2);
lmf = {'train', 'test', 'REt', 'T', 'REx', 'X', 'blank'};
trialtypes2probe = struct();
for lm = 1:numel(lmf)
    switch lmf{lm}
        case 'train'
            tempyu = traintrialtypes;
        case 'test'
            tempyu = traintrialtypes;
        case 'REt'
            tempyu = [1105, 1109];
        case 'T'
            tempyu = [105, 109];
        case 'REx'
            tempyu = [1201, 1299];
        case 'X'
            tempyu = [101, 1201, 1299];
        case 'blank'
            tempyu = 0;
    end
    trialtypes2probe.(lmf{lm}) = tempyu;
end
HR_SVMtrainICRC = struct();
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for lm = 1:numel(lmf)
            Ntt2probe = length(trialtypes2probe.(lmf{lm}));
            HR_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm}) = ...
                NaN(Ntt2probe,Ntt,Nsplits,Nsessions);
        end
    end
end

tic
if loadallsvm
    SVMtrainICRCagg = struct();
end
Nneuronsperarea = zeros(Nsessions, numel(visareas));
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    disp(whichvisarea)
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        disp(whichICblock)
        for ises = 1:Nsessions
            pathsvm = [pathsv nwbsessions{ises} filesep];
            svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
            if ~exist(svmfn, 'file')
                disp([svmfn ' does not exist'])
                continue
            end
            load(svmfn, 'SVMtrainICRC')
            if loadallsvm
            SVMtrainICRCagg.(whichICblock)(ises).(whichvisarea) = SVMtrainICRC;
            %SVMtrainICRCagg(ises).(whichICblock).(whichvisarea) = SVMtrainICRC;
            end
            Nneuronsperarea(ises,a) = nnz(SVMtrainICRC.neu2anal);
            for lm = 1:numel(lmf)
                tempyu = trialtypes2probe.(lmf{lm});
                switch lmf{lm}
                    case 'train'
                        tempyind = SVMtrainICRC.spkcnt.traintrialinds;
                        tempy = SVMtrainICRC.trialorder(tempyind);
                        tempysvm = SVMtrainICRC.spkcnt.train.label;
                    case 'test'
                        tempyind = SVMtrainICRC.spkcnt.testtrialinds;
                        tempy = SVMtrainICRC.trialorder(tempyind);
                        tempysvm = SVMtrainICRC.spkcnt.test.label;
                    otherwise
                        tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                        tempy = repmat( reshape(SVMtrainICRC.trialorder(tempyind),[],1),1,Nsplits) ;
                        tempysvm = SVMtrainICRC.spkcnt.all.label(tempyind,:);
                end
                if ~isequal(size(tempy), size(tempysvm))
                    error('check code')
                end
                for iy = 1:numel(tempyu)
                    trialsoy = tempy==tempyu(iy);
                    for jy = 1:numel(traintrialtypes)
                        trialsinclass = tempysvm==traintrialtypes(jy);
                        HR_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm})(iy,jy,:,ises) = ...
                            sum(trialsoy & trialsinclass, 1)./sum(trialsoy,1);
                    end
                end
            end

        end % ises
        toc
    end % b
end % a
toc

save([pathsv 'HR_SVMtrainICRC_' preproc '_selectareas_agg.mat'], 'SVMtrainICRCagg', 'Nneuronsperarea', 'HR_SVMtrainICRC', '-v7.3')
save(['G:\My Drive\RESEARCH\ICexpts_revision23\openscope_HR_SVMtrainICRC_' preproc '_agg.mat'], 'SVMtrainICRCagg', 'Nneuronsperarea', 'HR_SVMtrainICRC', '-v7.3')

%% compare areas
fs = 10;
discardbelowNneurons = 50;
whichICblock = 'ICwcfg1';
xtl = {'I_C_1', 'L_C_1', 'L_C_2', 'I_C_2'};
ytl = {'T_R_E_1', 'T_R_E_2'};

figure('Position',[100 100 800 400])
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, length(traintrialtypes)^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(length(traintrialtypes))),:), 1);
    sesoi = Nneuronsperarea(:,a)>=discardbelowNneurons;
    % sesoi = hrvec>1/length(traintrialtypes);
    % sesoi = true(size(tempHR,3),1);
    p = signrank(hrvec(sesoi) - 1/length(traintrialtypes));
    subplot(2,3,a)
    imagesc(nanmean(tempHR(:,:,sesoi), 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(xtl), 'YTickLabel', xtl)
    title(sprintf('%s %s %.2f\np=%.4f', whichvisarea, neuopt, mean(hrvec(sesoi)), p) )
end
colormap redblue

figure('Position',[100 500 800 400])
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, length(traintrialtypes)^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(length(traintrialtypes))),:), 1);
    sesoi = Nneuronsperarea(:,a)>=discardbelowNneurons;
    % sesoi = hrvec>1/length(traintrialtypes);
    % sesoi = true(size(tempHR,3),1);
    subplot(2,3,a)
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt, 3 ));
    imagesc(nanmean(tempHR(:,:,sesoi), 3))
    caxis([0 0.5]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(ytl), 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    p = signrank(infscore(sesoi));
    title(sprintf('%s %s IC-RC %.2f\np=%.4f', whichvisarea, neuopt, mean(infscore(sesoi)), p) )
end
colormap redblue


figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, length(traintrialtypes)^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(length(traintrialtypes))),:), 1);
    subplot(2,6,a)
    plot(Nneuronsperarea(:,a), hrvec, 'o')
    xlabel('# Neurons')
    ylabel('Cross-Validation Accuracy')
    title(sprintf('%s %.2f\np=%.4f', whichvisarea, mean(hrvec), p) )

    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt, 3 ));
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    subplot(2,6,6+a)
    plot(Nneuronsperarea(:,a), infscore, 'o')
    xlabel('# Neurons')
    ylabel('Inference Score')
    title(sprintf('%s IC-RC %.2f\np=%.4f', whichvisarea, mean(infscore), p) )
end
