maxratio = 0.50;
maxratio = 0.25;
maxratio = input('Enter which limiting factor you wan to look at (0.25 0.50 [default] 0.75 1.00) : ');
if length(maxratio) == 0
  maxratio = 0.50;
end
if length(intersect(maxratio,[0.25 0.50 0.75 1.00])) ~= 1
  error('can only use maxratio of 0.25, 0.50, 0.75, 1.00')
end
fprintf(1,'deltaT/deltaWV/deltaO3/deltaSKT updated according to uncertainty : max delta factor allowed = %8.3f \n',maxratio);
