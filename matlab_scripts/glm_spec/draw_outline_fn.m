function draw_outline_fn(time_ax, freq_ax, mask)
% Draw outline around the True regions of mask, on the current axes.
% time_ax: 1×nT vector of time bin centers
% freq_ax: 1×nF vector of frequency bin centers
% mask:    nF×nT logical (frequency × time, matching imagesc convention)

if ~any(mask(:))
    return
end

% Half bin widths (assume regular spacing; use first interval if not)
dt = mean(diff(time_ax))/2;
df = mean(diff(freq_ax))/2;

% Pad mask with zeros around the edges so we catch boundary transitions
m = false(size(mask) + 2);
m(2:end-1, 2:end-1) = mask;

% For each True cell in original mask, check its 4 neighbours.
% If a neighbour is False (or off-edge), draw the corresponding edge.
[nF, nT] = size(mask);
for fi = 1:nF
    for ti = 1:nT
        if ~mask(fi, ti)
            continue
        end
        tc = time_ax(ti);
        fc = freq_ax(fi);

        % Top edge (towards higher freq) - neighbour at fi+1
        if ~m(fi+2, ti+1)
            plot([tc-dt, tc+dt], [fc+df, fc+df], 'k', 'LineWidth', 1.5)
        end
        % Bottom edge (towards lower freq) - neighbour at fi-1
        if ~m(fi, ti+1)
            plot([tc-dt, tc+dt], [fc-df, fc-df], 'k', 'LineWidth', 1.5)
        end
        % Right edge (towards later time) - neighbour at ti+1
        if ~m(fi+1, ti+2)
            plot([tc+dt, tc+dt], [fc-df, fc+df], 'k', 'LineWidth', 1.5)
        end
        % Left edge (towards earlier time) - neighbour at ti-1
        if ~m(fi+1, ti)
            plot([tc-dt, tc+dt*0], [fc-df, fc+df], 'k', 'LineWidth', 1.5)
        end
    end
end
end