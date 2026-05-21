function fig = plot_eigs_unit_circle(lambda_vec, fig)
%==========================================================================
% Plots eigenvalues in the complex plane together with the unit circle.
%
% INPUTS:
%
%   lambda_vec     Eigenvalue vector                       [-]
%   fig            Figure handle (optional)               [-]
%
% OUTPUT:
%
%   fig            Figure handle                          [-]
%
%==========================================================================

    if nargin < 2 || isempty(fig)
        fig = figure('Position', [100, 100, 1000, 800]);
    else
        figure(fig);
    end

    hold on;
    grid on;
    box on;

    %---------------- Unit circle ----------------%
    theta_vec = linspace(0, 2*pi, 1000);

    plot(cos(theta_vec), sin(theta_vec), ...
        '--', ...
        'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.5);

    %---------------- Axes ----------------%
    plot([-1.2 1.2], [0 0], '-', ...
        'Color', [0.7 0.7 0.7], ...
        'LineWidth', 1);

    plot([0 0], [-1.2 1.2], '-', ...
        'Color', [0.7 0.7 0.7], ...
        'LineWidth', 1);

    %---------------- Eigenvalues ----------------%
    n_eigs = length(lambda_vec);

    % Same color for conjugate pairs
    pair_colors = lines(ceil(n_eigs/2));

    h = gobjects(n_eigs,1);
    labels = cell(n_eigs,1);

    for k = 1:n_eigs

        a = real(lambda_vec(k));
        b = imag(lambda_vec(k));

        %---------------- Clean numerical noise ----------------%
        tol = 1e-12;

        if abs(a) < tol
            a = 0;
        end

        if abs(b) < tol
            b = 0;
        end

        %---------------- Angle ----------------%
        theta_k = atan2(b, a);
        theta_deg = rad2deg(theta_k);

        if abs(theta_deg) < 1e-10
            theta_deg = 0;
        end

        %---------------- Plot point ----------------%
        color_idx = ceil(k/2);

        h(k) = plot(a, b, '.', ...
            'Color', pair_colors(color_idx,:), ...
            'MarkerSize', 30);

        %---------------- Scientific notation strings ----------------%
        a_str     = sci_notation_latex(a);
        b_abs_str = sci_notation_latex(abs(b));
        theta_str = sci_notation_latex(theta_deg);

        % Imaginary sign
        if b >= 0
            imag_sign = '+';
        else
            imag_sign = '-';
        end

        %---------------- Legend label ----------------%
        labels{k} = sprintf([ ...
            '$\\lambda_{%d} = %s %s %si$,  ' ...
            '$\\theta_{%d} = %s\\,\\mathrm{deg}$'], ...
            k, ...
            a_str, ...
            imag_sign, ...
            b_abs_str, ...
            k, ...
            theta_str);

    end

    %---------------- Legend ----------------%
    lgd = legend(h, labels, ...
        'Interpreter', 'latex', ...
        'FontSize', 20, ...
        'Location', 'eastoutside');

    lgd.Box = 'off';

    % Add spacing between legend entries
    lgd.ItemTokenSize = [25, 18];

    %---------------- Formatting ----------------%
    axis equal;

    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);

    xlabel('$\mathrm{Re}(\lambda)$', ...
        'Interpreter', 'latex', ...
        'FontSize', 20);

    ylabel('$\mathrm{Im}(\lambda)$', ...
        'Interpreter', 'latex', ...
        'FontSize', 20);

    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 1.5;

end

%==========================================================================
% Helper function:
% Converts numbers into LaTeX-friendly scientific notation
%==========================================================================

function str = sci_notation_latex(x)

    %---------------- Clean tiny values ----------------%
    if abs(x) < 1e-12
        x = 0;
    end

    %---------------- Zero case ----------------%
    if x == 0
        str = '0';
        return;
    end

    %---------------- Scientific notation ----------------%
    exponent = floor(log10(abs(x)));
    mantissa = x / 10^exponent;

    % Avoid scientific notation for moderate numbers
    if abs(exponent) <= 1

        str = sprintf('%.2f', x);

    else

        str = sprintf('%.2f \\times 10^{%d}', ...
            mantissa, exponent);

    end

end