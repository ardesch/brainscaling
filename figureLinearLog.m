function f = figureLinearLog(y, x, slope, intercept, label_y, label_x, varargin)
% This function creates a scatter plot with logarithmic x axis and regular
% y axis, with the option to add custom regression lines and confidence
% bands.

    % Defaults
    extra_line = [];
    mdl = [];
    
    % Parse optional arguments
    while ~isempty(varargin)
        if numel(varargin) == 1
            error('lscatter:missing_option', ...
                'Optional arguments must come in pairs.');
        end

        switch lower(varargin{1})
            case 'addline'
                assert(isnumeric(varargin{2}), 'addLine must be numeric');
                extra_line = varargin{2};
            case 'confidencebands'
                mdl = varargin{2}; % should be linear model with log x
        end

        varargin(1:2) = [];

    end

    % Plot x and y
    f = figure('Color', 'white', 'position', [200,200,300,300], 'Renderer', 'painters');
    s = scatter(x, y, 75, 'k', 'filled', 'MarkerFaceAlpha', 0.8);
    xlabel(label_x);
    ylabel(label_y);
    
    % Add regression line
    lims = [find(x == min(x)), find(x == max(x))];
    line_y = log(x(lims) .^ slope) + intercept; % linear-log plot: y = a*log(x)+b, so y = log(x^a) + b
    l = line(x(lims), line_y);
    l.LineStyle = '--';
    l.Color = 'k';
    
    % Add optional line
    if ~isempty(extra_line)
        
        hold on;
        tmp = log(x(lims) .^ extra_line) + intercept; % use same intercept as main line
        delta = line_y(1) - tmp(1);
        line2_y = tmp + delta;
        l2 = line(x(lims), line2_y);
        l2.LineStyle = '--';
        l2.Color = [0.5, 0.5, 0.5];        
        
    end
    
    % Add optional confidence bands
    if ~isempty(mdl)
        
        hold on;
        logxmin = min(log(x));
        logxmax = max(log(x));
        logxrange = logxmax - logxmin;
        xlogci = [logxmin:logxrange/1000:logxmax]'; % define some points for x axis on log scale
        [~, yci] = predict(mdl, xlogci);
        xci = exp(xlogci);
        patchx = [xci; xci(end:-1:1)];
        patchy = [yci(:,1); yci(end:-1:1,2)];
        patch('XData', patchx, 'YData', patchy, ...
            'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
    
    end
    
    % Put the data points back on top
    uistack(s,'top')
    
    % Set axis to log scale
    set(gca, 'XScale', 'log')
   
    % Set axis limits
    xlim([min(x).*0.75, max(x).*1.25]);
    y_spacing = 0.10 * min(y);
    ylim([min(y) - y_spacing, max(y) + y_spacing]);
    
    % Set font
    set(gca, 'FontSize', 14);
    set(gca, 'FontName', 'Helvetica');    

end