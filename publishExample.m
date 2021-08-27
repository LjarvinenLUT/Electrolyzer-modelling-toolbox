%% A script for publishing example.m in both HTML and latex

% HTML
options = struct('format','html','outputDir',[pwd '\Example']);
publish('example.m',options)

% % LaTeX
% options.format = 'latex';
% options.stylesheet = 'publish2latex.xsl';
% publish('example.m',options)
