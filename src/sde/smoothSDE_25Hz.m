% Pre-processing data from .mat PRH files for smoothSDE analysis for Atl
% BRS
%
% Note: this script requires animaltags tool kit, available at:
% https://github.com/animaltags/tagtools_matlab

% tag IDs to process
tagz = [
    'zc17_234a';
    'zc19_218a';
    'zc20_232a';
    'zc22_219a' % this one is real MFAS, rest are all scaled
    ];


% settings for dive detection
mindepth = 50;
at_surf = 2;


% loop over tags and process one by one, saving results in a table()
all_whales = table(); % empty place holder
for t = 1:size(tagz,1)
    clear A Aw fs head M Mw p pitch roll D
    % read in data from PRH file (if you have 12.5 and 25 Hz USE 25)
    load(['data/prh/', tagz(t,:), 'prh.mat']);
    % dive detections: using settings @ top of script, findall = 0
    D = find_dives(p, fs, mindepth, at_surf);
    mag_jerk = njerk(M, fs);
    msa1 = msa(A, 1);
   
    % convert  heading to degrees
    head_deg = head .* 180/pi;

    % time stamps in HOURs since start of record
    % place in mid-point of the 3-second time window
    time_hr = [1:size(p, 1)]' / fs / 3600;

    % whale ID
    whaleID = repmat(tagz(t,:), size(time_hr, 1), 1);

    % save results in a table
    this_tab = table(whaleID, time_hr, p, ...
        pitch, roll, head_deg, ...
        msa1, mag_jerk);

    % combine all data so far
    all_whales = vertcat(all_whales, this_tab);
end

% write results to text file
writetable(all_whales, 'data/AtlBRS_Zc_smoothSDE_data_25Hz.csv');

