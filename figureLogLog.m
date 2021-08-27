function f = figureLogLog(y, x, slope, intercept, label_y, label_x, varargin)
% This function creates a scatter plot with logarithmic x axis and y axes, 
% with the option to add custom regression lines and confidence bands.

    % Defaults
    extra_line = [];
    mdl = [];
    datapoints = [];
    
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
            case 'plotbustomband'
                datapoints = varargin{2};
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
    line_y = exp(intercept) .* x(lims) .^ slope;
    l = line(x(lims), line_y);
    l.LineStyle = '--';
    l.Color = 'k';
    
    % Add optional line
    if ~isempty(extra_line)
        
        hold on;
        tmp = exp(intercept) .* x(lims) .^ extra_line; % use same intercept as main line
        delta = log(line_y(1)) - log(tmp(1));
        
        % scale such that additional line goes 45 degrees
        % separate cases for when extra line slope is higher or lower
        % (positive or negative allometric scaling)
        if extra_line < slope % positive scaling
            
            % find x for the extra_line where y = max(line_y)
            line2_x_min = min(x);
            line2_y_min = exp(intercept+delta) .* min(x) .^ extra_line;
            line2_x_max = (max(y)/exp(intercept+delta)) .^ (1/extra_line);
            line2_y_max = exp(intercept+delta) .* line2_x_max .^ extra_line;
            
            axis_lims_x = [min(x) .* 0.75, line2_x_max .* 1.25];
            axis_lims_y = [min(y) .* 0.75, line2_y_max .* 1.25];
            
        else % negative or isometric scaling
            
            line2_y = exp(intercept+delta) .* x(lims) .^ extra_line; % add difference in raw data y to intercept
            line2_y_min = min(line2_y);
            line2_y_max = max(line2_y);
            axis_lims_x = [min(x) .* 0.75, max(x) .* 1.25];
            axis_lims_y = [min(y) .* 0.75, max(line2_y)];
            line2_x_min = min(x);
            line2_x_max = max(x);
    
        end
        
        l2 = line([line2_x_min, line2_x_max], [line2_y_min, line2_y_max]);
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
        [~, ylogci] = predict(mdl, xlogci);
        yci = exp(ylogci);
        xci = exp(xlogci);
        patchx = [xci; xci(end:-1:1)];
        patchy = [yci(:,1); yci(end:-1:1,2)];
        patch('XData', patchx, 'YData', patchy, ...
            'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);    
        
    end
    
    % Add optional confidence bands based on x and y values specified in
    % datapoints
    if ~isempty(datapoints)
        
        hold on;
        xci = datapoints(:,1);
        yci = datapoints(:,2:3);
        patchx = [xci; xci(end:-1:1)];
        patchy = [yci(:,1); yci(end:-1:1,2)];
        patch('XData', patchx, 'YData', patchy, ...
            'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);    
        
    end
    
    % Put the data points back on top
    uistack(s,'top')
    
    % Set axis to log scale
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
   
    % Set axis limits
    xlim(axis_lims_x);
    ylim(axis_lims_y);

    % Set font
    set(gca, 'FontSize', 14);
    set(gca, 'FontName', 'Helvetica');    
    
end