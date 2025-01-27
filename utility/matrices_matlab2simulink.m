model = 'deltino';

fprintf('Starting creating Simulink block for matrix M...\n')
tic; 
matlabFunctionBlock([model '/M'], M); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

fprintf('Starting creating Simulink block for matrix C...\n')
tic; 
matlabFunctionBlock([model '/C'], C); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

fprintf('Starting creating Simulink block for vector G...\n')
tic; 
matlabFunctionBlock([model '/G'], G); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

fprintf('Starting creating Simulink block for matrix A...\n')
tic; 
matlabFunctionBlock([model '/A'], A); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

fprintf('Starting creating Simulink block for matrix S...\n')
tic; 
matlabFunctionBlock([model '/S'], S); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

fprintf('Starting creating Simulink block for matrix S_dot...\n')
tic; 
matlabFunctionBlock([model '/S_dot'], S_dot); 
duration = toc;
fprintf('... duration = %d s\n\n', duration)

