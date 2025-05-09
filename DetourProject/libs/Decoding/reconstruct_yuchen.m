function [pdf,tbins,nact] = reconstruct_yuchen(ts, te, tc, cells, varargin)
    
%     ts: start time
%     te: stop time
%     tc: place field. spatial bins X no.of cells
%     cells : cells(i).time
%     nact is number of active units in each time bin
% % %   example to modify parameter:
% % %   reconstruct(ts, te, tc, cells, 'tau', 0.02)
% % %   now args.tau is overwrited to 0.02...

 
    args.tau = 0.25;
    args.percent_overlap =0; % originally was 0
    args.max_pdf_comp_size = 5000; % originally was 5000, then 2000
 
    
    args = parseArgs(varargin, args);
    tbins = get_tbins(ts, te, args.tau, args.percent_overlap);
    spike_counts= zeros(length(cells), size(tbins,1));
    for i=1:length(cells)
        spike_counts(i,:) = histcounts_row_kf(cells(i).time, tbins);
    end
    
    
    % get number of active units in each time bin
    acttbin = spike_counts >= 1;
    nact = sum(acttbin,1);

    n_dir = size(tc,3);
    pdf = nan(size(tc,1), size(tbins,1), n_dir);
    for dir = 1:n_dir
        if size(spike_counts,2)>args.max_pdf_comp_size
            n_section = ceil(size(spike_counts,2)/args.max_pdf_comp_size);
            disp(['Requested PDF is too big, cutting into ', num2str(n_section), ' sections']);
            for i=0:n_section-1
                istart = i*args.max_pdf_comp_size+1;
                iend = min([(i+1)*args.max_pdf_comp_size, size(spike_counts,2)]);
%                 size(tc)
%                 size(spike_counts)
                pdf_short = parameter_estimation_simple(args.tau, tc(:,:,dir), spike_counts(:,istart:iend));

                switch i
                    case 0
                        pdf_temp = pdf_short;
                    otherwise
                        pdf_temp = [pdf_temp pdf_short]; %#ok
                end
            end
        else
            pdf_temp = parameter_estimation_simple(args.tau, tc(:,:,dir), spike_counts);
        end
        pdf = pdf_temp;
    end
    
end



function tbins = get_tbins(ts, te, tau, percent_overlap)
tbins = ts:tau*(1-percent_overlap):(te-tau);
if tbins(end) + tau < te
   tbins = [tbins,tbins(end)+tau*(1-percent_overlap)]; 
end

if numel(tbins)==1
    tbins = [tbins, tbins+tau];
end
tbins = [tbins', tbins'+tau];
end

