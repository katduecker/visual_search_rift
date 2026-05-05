function draw_outline_fn(time_ax, freq_ax, mask)
% Draw outlines around contiguous True regions in mask.
% time_ax: 1×nT bin centers
% freq_ax: 1×nF bin centers
% mask:    nF×nT logical (rows = freq, cols = time, matching imagesc convention after axis xy)

if ~any(mask(:)), return; end

% Force row vectors
time_ax = time_ax(:)';
freq_ax = freq_ax(:)';

% Half bin widths from the actual axes
dt = mean(diff(time_ax)) / 2;
df = mean(diff(freq_ax)) / 2;

[nF, nT] = size(mask);

for fi = 1:nF
    for ti = 1:nT
        if ~mask(fi, ti), continue; end

        tL = time_ax(ti) - dt;   % left edge of this bin in time
        tR = time_ax(ti) + dt;   % right edge
        fB = freq_ax(fi) - df;   % bottom edge in freq
        fT = freq_ax(fi) + df;   % top edge

        % Bottom edge: draw if the cell below (lower freq) is outside cluster or off-grid
        if fi == 1 || ~mask(fi-1, ti)
            line([tL, tR], [fB, fB], 'Color', 'k', 'LineWidth', 1.5)
        end
        % Top edge
        if fi == nF || ~mask(fi+1, ti)
            line([tL, tR], [fT, fT], 'Color', 'k', 'LineWidth', 1.5)
        end
        % Left edge
        if ti == 1 || ~mask(fi, ti-1)
            line([tL, tL], [fB, fT], 'Color', 'k', 'LineWidth', 1.5)
        end
        % Right edge
        if ti == nT || ~mask(fi, ti+1)
            line([tR, tR], [fB, fT], 'Color', 'k', 'LineWidth', 1.5)
        end
    end
end
end