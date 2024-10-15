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

% settings for decimation of data 
ave_win = 3; % output will be one sample per ave_win seconds
% settings for dive detection
mindepth = 50;
at_surf = 2;


% loop over tags and process one by one, saving results in a table()
all_whales = table(); % empty place holder
for t = 1:size(tagz,1)
    clear A Aw fs head M Mw p pitch roll D
    % read in data from PRH file (if you have 12.5 and 25 Hz USE 25)
    % (why: decimation can only be by integer factor and we want 3 Hz)
    load(['data/prh/', tagz(t,:), 'prh.mat']);
    % dive detections: using settings @ top of script, findall = 0
    D = find_dives(p, fs, mindepth, at_surf);
    % decimate data to 1 sample per 3 seconds, or fs = 1/3 Hz
    pitch_lo = decdc(pitch, ave_win * fs);
    head_lo = decdc(head, ave_win * fs);
    roll_lo = decdc(roll, ave_win * fs);
    p_lo = decdc(p, ave_win * fs);
    mag_jerk_mean = block_mean(njerk(M, fs), ave_win * fs);
    mag_jerk_lo = decdc(njerk(M, fs), ave_win * fs);
    msa_lo = msa(A, 1);
    msa_lo = decdc(msa_lo, ave_win * fs) ;
    
    % convert  heading to degrees
    head_lo_deg = head_lo .* 180/pi;

    % time stamps in HOURs since start of record
    % place in mid-point of the 3-second time window
    time_hr = [1:size(p_lo, 1)]' * ave_win / 60 / 60 + (ave_win / 2 / 60 / 60);

    % which dive number?
    ID = ones(size(time_hr));
    % proportion of total time in dive
    diveprop = NaN * ID;
    max_depth = diveprop;
    for r = 1:size(time_hr, 1)
        if (time_hr(r) < (D.start(1) /60 /60)) || (time_hr(r) > (D.end(end) / 60 /60))
            ID(r) = NaN;
        else
            ID(r) = find(time_hr(r) > (D.start / 60 / 60), 1, 'last');
            diveprop(r) = (time_hr(r) - (D.start(ID(r)) / 60 / 60)) / ...
                ((D.end(ID(r)) - D.start(ID(r))) / 60 / 60);
            max_depth(r) = D.max(ID(r));
        end
    end

    % whale ID
    whaleID = repmat(tagz(t,:), size(time_hr, 1), 1);

    % save results in a table
    this_tab = table(whaleID, ID, time_hr, diveprop, p_lo, ...
        pitch_lo, roll_lo, head_lo_deg, max_depth, ...
        msa_lo, mag_jerk_lo, mag_jerk_mean);

    % combine all data so far
    all_whales = vertcat(all_whales, this_tab);
end

% write results to text file
writetable(all_whales, 'data/AtlBRS_Zc_smoothSDE_data.csv');

