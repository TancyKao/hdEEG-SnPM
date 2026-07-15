function [gpEEG] = func_mergeKDTstages(curEEG)

eyeopen_idx = find(contains(curEEG.stages, 'eyesopen'));
eyeclose_idx = find(contains(curEEG.stages, 'eyesclosed'));

eyeo_dat = mean(curEEG.data(:,eyeopen_idx,:,:),2,'omitnan');
eyec_dat = mean(curEEG.data(:,eyeclose_idx,:,:),2,'omitnan');

curEEG.data = cat(2,eyeo_dat, eyec_dat);
curEEG.stages = {'eyesopen','eyesclosed'};

gpEEG = curEEG;





