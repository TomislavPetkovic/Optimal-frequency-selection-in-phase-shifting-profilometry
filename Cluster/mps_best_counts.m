%% Selection of best fringe counts
% Use this script to find best possible fringe counts using exhaustive
% search.

clear all;
close all;


%% ADD PATHS
% Add paths to directories containing required functions.

%phase_shift_functions = 'c:\Users\Tomislav\Documents\3DTS\Matlab\Phase Shift';
%mps_functions = 'c:\Users\Tomislav\Documents\3DTS\Matlab\Multiple Phase Shift';
%optimal_mps_functions = 'C:\Users\Tomislav\Documents\3D-CODING\Matlab\Optimal MPS Frequencies';

phase_shift_functions = '/home/tpetkovic/Phase Shift';
mps_functions = '/home/tpetkovic/Multiple Phase Shift';
optimal_mps_functions = '/home/tpetkovic/Optimal MPS Frequencies';

s = pathsep;
p = [s path s];
if ~contains(p, [s phase_shift_functions s], 'IgnoreCase', true); path(phase_shift_functions, path); end
if ~contains(p, [s mps_functions s], 'IgnoreCase', true); path(mps_functions, path); end
if ~contains(p, [s optimal_mps_functions s], 'IgnoreCase', true); path(optimal_mps_functions, path); end


%% SET PARAMETERS
% We have to run exhaustive search for multiple fringe counts. Here we set
% parameters for each individual search.

%N_all =         {  2,   3,   4,   5,   6,   2,   3,   4,   5,   6,  2,  3,  4,  5,  6,   2,   3,   4,   5,   6};
%N_all =         {  2,   2,   2,   2,   2,   2,   2,   2,   2,   2,  2,  2,  2,  2,  2,   2,   2,   2,   2,   2};
%count_min_all = { 64,  64,  64,  64,  64,  80,  80,  80,  80,  80, 40, 40, 40, 40, 40,  50,  50,  50,  50,  50};
%count_max_all = {128, 128, 128, 128, 128, 160, 160, 160, 160, 160, 80, 80, 80, 80, 80, 100, 100, 100, 100, 100};

%N_all =         { 5,   5,   5,  6,   6,   6};
%count_min_all = {40,  50,  80, 40,  50,  80};
%count_max_all = {80, 100, 160, 80, 100, 160};

%N_all =         {  2,   3,   4,   2,   3,   4,  2,  3,  4,   2,   3,   4};
%count_min_all = { 64,  64,  64,  80,  80,  80, 40, 40, 40,  50,  50,  50};
%count_max_all = {128, 128, 128, 160, 160, 160, 80, 80, 80, 100, 100, 100};

%N_all =         {  2,   2,  2,   2};
%count_min_all = { 64,  80, 40,  50};
%count_max_all = {128, 160, 80, 100};

%N_all =         {  5,   5,  5,   5};
%count_min_all = { 64,  80, 40,  50};
%count_max_all = {128, 160, 80, 100};

N_all =         {  6};
count_min_all = { 80};
count_max_all = {160};

sz = size(N_all);

d_best = cell(sz);
count_best = cell(sz);
lambda_best = cell(sz);
width_best = cell(sz);

d_worst = cell(sz);
count_worst = cell(sz);
lambda_worst = cell(sz);
width_worst = cell(sz);

num_threads = 20;

basename = 'mps_best_counts';

diaryname = [basename '_diary.txt'];
diary(diaryname);
diary ON;


%% START PARALLEL POOL
% Start parallel pool. This function should be run once.

if isempty(gcp('nocreate'))
    myPool = parpool(num_threads);
else
    myPool = gcp;
end


%% RUN SEARCH
% This is the outer loop for exhaustive search.

jstart = 0;
if 0 < jstart
    fprintf('Loading data for slot %d.\n', jstart);
    load([basename '_' num2str(jstart) '.mat'], 'd_best', 'count_best', 'lambda_best', 'width_best', 'd_worst', 'count_worst', 'lambda_worst', 'width_worst');
    
    assert( all(size(d_best) == sz) );
    assert( all(size(count_best) == sz) );
    assert( all(size(lambda_best) == sz) );
    assert( all(size(width_best) == sz) );
    
    assert( all(size(d_worst) == sz) );
    assert( all(size(count_worst) == sz) );
    assert( all(size(lambda_worst) == sz) );
    assert( all(size(width_worst) == sz) );
end
jstart = jstart + 1;

for j = jstart : numel(N_all)
    
    %% SET PARAMETERS
    % Exaustive search simply generates all counts within the specified
    % interval. Then for each count a minimal constellation distance is
    % computed. The best fringe counts are ones which produce the largest
    % constellation distance.
    % Regarding parameters note that all valid combinations will be generated
    % and tested which may take a lot of memory and be quite slow for large N.
    
    N = N_all{j}; % Number of fringe counts (number of frequencies).
    assert(2 <= N);
    
    count_min = count_min_all{j}; % Minimal number of fringe counts per screen.
    count_max = count_max_all{j}; % Maximal number of fringe counts per screen.
    assert( (0 < count_min) && (0 < count_max) && (count_min <= count_max) );
    
    fprintf('\n\n------ NEW EXHAUSTIVE SEARCH ------\n');
    fprintf('Number of fringes: %d\n', N);
    fprintf('Minimal fringe count per screen: %d\n', count_min);
    fprintf('Maximal fringe count per screen: %d\n', count_max);
    
    filename_j = [basename '_N' num2str(N) '_min' num2str(count_min) '_max' num2str(count_max) '.mat'];
    
    
    %% GENERATE ALL VALID COUNTS
    % First we generate all N tuples which are valid counts. An N tuple is
    % a valid count if all of its elements are unique and if the greatest
    % common divisor of its elements is 1.
    
    % First generate all counts and compute the upper limit on the number
    % of combinations.
    V = count_min : count_max;
    V_len = length(V);
    num_tuples = nchoosek(V_len, N);
    
    mem_available = NaN;
    if ~ispc
        [status, result] = system('free -b | grep Mem');
        if 0 == status
            mem = str2double(regexp(result, '[0-9]*', 'match'));
            mem_available = round(mem(1) * 0.8);
        end
    else
        mem = memory;
        mem_available = mem.MemAvailableAllArrays;
    end
    if ~isnan(mem_available)
        mem_required = N * num_tuples * 8;
        fprintf('Need %d bytes to store all combinations.\n', mem_required);
        if mem_required > mem_available
            warning('Generating all combinations may use swap memory!');
            pause;
        end
    end
    
    % Generate all combinations.
    max_abs_V = max(abs(V(:)));
    if 2^15 - 1 > max_abs_V
        V = int16(V);
    elseif 2^31 - 1 > max_abs_V
        V = int32(V);
    end    
    counts_all = nchoosek(V, N);
    counts_all = counts_all.';
    M = size(counts_all, 2);
    fprintf('Generated %d combinations.\n', M);
    
    % Then retain only valid combinations; we us in-place processing.
    K = 0;
    for i = K + 1 : M
        count = counts_all(:, i);
        assert( numel(unique(count)) == N );
        if 1 == mps_gcd(count)
            K = K + 1;
            assert(K <= i);
            counts_all(:, K) = count;
        end
    end
    counts_all = counts_all(:, 1:K);
    fprintf('Of %d combinations only %d are valid.\n', M, K);
    
    
    %% EVALUATE EACH VALID COUNT
    % For each fringe count generate the correspoding periods (wavelengths)
    % and compute the constellation half-distance. We do this using a
    % parallel pool of workers to speed up the search.
        
    % Restart parallel pool if needed. 
    if isempty(gcp('nocreate'))
        myPool = parpool(num_threads);
    else
        myPool = gcp;
    end
    
    % Strat time.
    tstart = tic;
    
    % Distribute work to parallel workers.
    q = parallel.pool.PollableDataQueue;
    nw = myPool.NumWorkers;
    assert( nw == num_threads );
    F(1:nw) = parallel.FevalFuture;    
    step = ceil(K / nw);
    i2_prev = 0;
    total = 0;
    for i = 1 : nw
        i1 = (i-1) * step + 1;
        assert( i1 - 1 == i2_prev );
        
        i2 = i1 + step - 1;
        if (i2 > K); i2 = K; end
        
        fprintf('Sending task %d to workers.\n', i);
        
        F(i) = parfeval(@mps_best_counts_min_max, 2, counts_all(:, i1:i2), q);
        
        total = total + (i2 - i1 + 1);
        i2_prev = i2;
    end
    assert( total == K );
    
    fprintf('Distributed all tasks to parallel workers.\n');
    
    % Wait for all workers to complete and then combine outputs.
    out_min = {};
    out_max = {};
    total_time = 0;
    num_processed = 0;
    timeout = 60 * 5; % Wait timeout in seconds.    
    while ~all([F.Read])
        % Try to fetch results until timeout.
        [idx, tmp_min, tmp_max] = fetchNext(F, timeout);        
        if ~isempty(idx)            
            if isempty(out_min)
                out_min = tmp_min;
                out_max = tmp_max;
            else
                [out_min, out_max] = mps_best_counts_combine(out_min, out_max, tmp_min, tmp_max);
            end
            duration = seconds( F(idx).FinishDateTime - F(idx).StartDateTime );
            total_time = total_time + duration;
            
            fprintf('Fetched results from task %d.\n', idx);
        end
        
        % Pool the queue for the number of processed items so far.
        have_messages = true;
        read_messages = false;
        while have_messages
            [value, have_messages] = poll(q, 0.01);
            if ~isempty(value)
                num_processed = num_processed + value;
                read_messages = true;
            end
        end
        if read_messages
            fprintf('Processed a total of %d items out of %d (%.1f %%).\n', num_processed, K, 100*num_processed/K);
            remaining = toc(tstart) * (K - num_processed) / num_processed;
            if 60 > remaining
                fprintf('Estimated time to completion %.1f seconds.\n', remaining);
            elseif 3600 > remaining
                remaining_min = floor(remaining/60);
                remaining_sec = remaining - remaining_min * 60;
                fprintf('Estimated time to completion %d minutes %.1f seconds.\n', remaining_min, remaining_sec);
            elseif 86400 > remaining
                remaining_hours = floor(remaining/3600);
                remaining_min = floor((remaining - remaining_hours * 3600)/60);
                remaining_sec = remaining - remaining_hours * 3600 - remaining_min * 60;
                fprintf('Estimated time to completion %d hours %d minutes %.1f seconds.\n', remaining_hours, remaining_min, remaining_sec);
            else
                remaining_days = remaining / 86400;
                fprintf('Estimated time to completion %.1f days.\n', remaining_days);
            end
        end
    end
    elapsed = toc(tstart);
    
    fprintf('Total running time of %d thread(s) is %.1f seconds.\n', nw, total_time);
    fprintf('Actual running time is %.1f seconds.\n', elapsed);
    
    % Assign final results.
    distance_max = out_max{1};
    lambda_max = out_max{2};
    width_max = out_max{3};
    counts_max = out_max{4};
    
    distance_min = out_min{1};
    lambda_min = out_min{2};
    width_min = out_min{3};
    counts_min = out_min{4};
    
    
    %% PRINT BEST AND WORST SOLUTIONS
    % Output the best solution according to the constellation half-distance
    % only. Note that other metric may be used to select the best solution. For
    % example higher counts have better properties with regard to
    % inter-reflections in the scene but always produce worse half-distance due
    % to constellation beeing more dense; therefore, a trade-off between having
    % higher fringe count per screen and having lower constellation
    % half-distance is possible.
    
    d_best{j} = distance_max;
    count_best{j} = counts_max;
    lambda_best{j} = lambda_max;
    width_best{j} = width_max;
    
    fprintf('\nThe best found solution for counts in [%d,%d] interval:\n  fringe counts %d', count_min, count_max, count_best{j}(1));
    for i = 2 : N; fprintf(':%d', count_best{j}(i)); end
    fprintf('\n  periods %d', lambda_best{j}(1));
    for i = 2 : N; fprintf(':%d', lambda_best{j}(i)); end
    fprintf('\n  width %d\n  minimal half-distance %.4f rad (%.2f deg)\n', width_best{j}, d_best{j}, d_best{j}*180/pi);
    
    d_worst{j} = distance_min;
    count_worst{j} = counts_min;
    lambda_worst{j} = lambda_min;
    width_worst{j} = width_min;
    
    fprintf('\nThe worst found solution for counts in [%d,%d] interval:\n  fringe counts %d', count_min, count_max, count_worst{j}(1));
    for i = 2 : N; fprintf(':%d', count_worst{j}(i)); end
    fprintf('\n  periods %d', lambda_worst{j}(1));
    for i = 2 : N; fprintf(':%d', lambda_worst{j}(i)); end
    fprintf('\n  width %d\n  minimal half-distance %.4f rad (%.2f deg)\n', width_worst{j}, d_worst{j}, d_worst{j}*180/pi);
    
    
    save([basename '_' num2str(j) '.mat'], 'd_best', 'count_best', 'lambda_best', 'width_best', 'd_worst', 'count_worst', 'lambda_worst', 'width_worst');

end

disp(' ');
disp('SEARCH COMPLETED');
diary OFF;