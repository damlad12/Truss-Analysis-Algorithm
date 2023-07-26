S = vertcat(Sx,Sy);
active_load = max(L);
j = length(X);

Lengths = zeros((2*j)-3, 1);

numcols = size(C, 2);
Tcoeffx = [];
Tcoeffy = [];


for col = 1:numcols
   columnx = zeros(j, 1);
   columny = zeros (j, 1);
   column = C(:, col);
   first = 0 ;
   second = 0; 
   for i = 1:length(column)
    if (column(i) == 1 && first == 0 )
        first = i;
    elseif (column(i) == 1 && first ~= 0 )
        second = i;
    end 
   end
    x1 = X(first);
    y1 = Y(first);
    p1 = [x1,y1];
    x2 = X(second);
    y2 = Y(second);
    p2 = [x2,y2];
    columnx(first) = (x2 - x1)/norm(p2-p1)  ;
    columnx(second) = (x1 - x2)/norm(p2-p1) ;
    columny(first) = (y2 - y1)/norm(p2-p1) ;
    columny(second) = (y1 - y2)/norm(p2-p1)  ; 
    Lengths(col) = norm(p2-p1);

    Tcoeffx = [Tcoeffx columnx];
    Tcoeffy = [Tcoeffy columny];  
  

end

Tcoeff = vertcat(Tcoeffx, Tcoeffy);

A = horzcat(Tcoeff, S);
T = inv(A) * L;

R = (T/active_load);

Pcrit = zeros((2*j)-3, 1);

for k = 1 : length(Lengths)
   Pcrit(k) = -(4338 * ((Lengths(k))^(-2.125))) ;
end

for k = 1 : length(Lengths)
      disp(Pcrit(k));
end

Wfailure = zeros((2*j)-3, 1);

for z = 1 : (length(Lengths))
Wfailure(z) = -Pcrit(z)/R(z);
end


max_load_neg = max(Wfailure(Wfailure < 0));
max_load = abs(max(Wfailure(Wfailure < 0)));



failing_member = 0; 
for z = 1 : length(Lengths)
if (Wfailure(z) == max_load_neg)
    failing_member = z;
end
end
disp(max_load);
disp(failing_member);
disp(R(failing_member));


fprintf('%% EK301, Section A1, Group DSK: Alexandra S.,Damla D., 4/7/2023.\n');
fprintf('Load: %d oz\n',active_load);
fprintf('%% Member forces in oz\n');

num_Sx = 0 ;
numcols2 = size(Sx, 2);
for col = 1:numcols2
   column = Sx(:, col);
   for i = 1:length(column)
       if  (column(i) == 1)
           num_Sx =  num_Sx + 1;
       end 
   end 
end 

num_Sy = 0 ;
numcols3 = size(Sy, 2);
for col = 1:numcols3
   column = Sy(:, col);
   for i = 1:length(column)
       if  (column(i) == 1)
           num_Sy =  num_Sy + 1;
       end 
   end 
end

for g = 1:(length(T) - (num_Sy + num_Sx))
    if(T(g)<0)
    fprintf('m%d: %.3f (C)\n', g, abs(T(g)));
    else 
     fprintf('m%d: %.3f (T)\n', g, abs(T(g)));
    end
end

fprintf('%% Reaction forces in oz\n');

o = 0;
numcols2 = size(Sx, 2);
for col = 1:numcols2
   column = Sx(:, col);
   for i = 1:length(column)
       if  (column(i) == 1)
           o = o + 1;
             fprintf('Sx%d: %.3f \n', o, abs(T(2*j- 3 + o)));
       end 
   end 
end 

o = 0;
numcols3 = size(Sy, 2);
for col = 1:numcols3
   column = Sy(:, col);
   for i = 1:length(column)
       if  (column(i) == 1)
           o = o + 1;
              fprintf('Sy%d: %.3f \n', o, abs(T(2*j- 3 + num_Sx + o)));
       end 
   end 
end 

sum_lengths = 0;
for i = 1:size(Lengths)
    sum_lengths = sum_lengths + Lengths(i);
end


cost = (j * 10 + sum_lengths);
fprintf('Cost of truss: $%d\n', cost);

max_load_cost_ratio = max_load/cost ;
fprintf('Theoretical max load/cost ratio in oz/$: %.4f\n', max_load_cost_ratio);
