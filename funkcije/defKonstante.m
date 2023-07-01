function defKonstante()
assignin('base', 'a', exp(i*2*pi/3));
evalin('base', ['global a'])
a = exp(i*2*pi/3);

assignin('base', 'T', [1 1 1; 1 a a*a; 1 a*a a]);
evalin('base', ['global T']);


end

