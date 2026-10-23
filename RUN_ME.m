function RUN_ME
% Compute the numerical example and export the figures.
packageRoot = fileparts(mfilename('fullpath'));
originalPath = path;
originalFigureVisibility = get(groot, 'defaultFigureVisible');
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
figureCleanup = onCleanup(@() set(groot, 'defaultFigureVisible', originalFigureVisibility)); %#ok<NASGU>
addpath(fullfile(packageRoot, '03_code'));
fprintf('MATLAB version: %s\n', version);
run_bam('reproduced_matlab');
outputFolder = fullfile(packageRoot, 'outputs', 'reproduced_matlab');
fprintf('\nResults folder:\n%s\n', outputFolder);
fprintf('Derivative-error figure: fig5_derivative_errors.png\n');
fprintf('Archived paper data in 04_results have not been overwritten.\n');
end
