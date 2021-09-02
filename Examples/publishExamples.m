%% A script for publishing example.m in both HTML and latex

% HTML
options = struct('format','html');
%%
publish('example.m',options)
%%
publish('funcCreationExample.m',options)

%% LaTeX
% options.format = 'latex';
% options.stylesheet = 'publish2latex.xsl';
% publish('example.m',options)
