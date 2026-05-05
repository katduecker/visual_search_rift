function draw_outline_fn(time_ax, freq_ax, m_vert, m_horz)
for i = 1:size(m_vert,1)
    for m = 2:size(m_vert,2)
        if m_vert(i,m)
            plot(time_ax(m-1:m)+0.025, ones(1,2)*freq_ax(i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end
for i = 2:size(m_horz,1)
    for m = 2:size(m_horz,2)
        if m_horz(i,m)
            plot(time_ax(m-1)*ones(1,2)+0.025, freq_ax(i-1:i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end