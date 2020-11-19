function lab2
    % Values to theta2 and theta4 are  provided in the question
    %------------part a-------------
    % angles are provided in the lab description
    theta2  = [40; 45; 50; 55; 60]; 
    theta4  = [70; 76; 83; 91; 100];
    link_ratios = get_link_ratios(theta2, theta4);
    [a,b,c,d] = get_link_lengths(link_ratios);
    fprintf("Crank: %.4f mm\n Coupler: %.4f mm\n Follower: %.4f mm\n  Fixed: %.4f mm\n",a,b,c,d);
    %----------- part b --------------
    transmission_angles  = get_transmission_angles(a,b,c,d,40,60,1);
    input_angles = 40:1:60;
    figure;
    title("Transmission angles vs Input angles");
    plot(input_angles, transmission_angles);
    xlabel("Input angles");
    ylabel("Transmission angles");
   % commenting  on the quality of the transmission angles
            if (all(transmission_angles >= 40)) || (all(transmission_angles <= 140))
                fprintf("All the transmission angles guarantee a smooth rotation\n");
            else
                fprintf("Some transmission angles do not guarantee a smooth rotation\n")
            end
     %-----------part c ---------------
     % [a1, b1] = applyRegression(theta2, theta4);
      % output_angles = arrayfun(@(x) x*b1 + a1 , input_angles);
      %structural errors
        A = sind(input_angles);
        B = cosd(input_angles)-link_ratios(1);
        C = link_ratios(3) - link_ratios(2)*cosd(input_angles);

        tan_1 = (A + sqrt(A.^2 + B.^2 -C.^2))./(B+C);
        %disp(tan_1)
        tan_2 = (A - sqrt(A.^2 + B.^2 -C.^2))./(B+C);
        %disp(tan_2)
        out_1 = 2*atand(tan_1);
        out_2 = 2*atand(tan_2);
      structural_errors = get_structural_errors(input_angles, out_2, link_ratios);
      figure;
      title("Structural errors Vs Input angles");
      plot(input_angles, structural_errors);
      xlabel("Input angles");
      ylabel("Structural errors");
end

function link_ratios  = get_link_ratios(theta2, theta4)
   % using the freudensteins method, this method computes the link ratios
   
   % To form a complete and functional matrix the two 1D arrays  should be
   % of the same length. If not, this raises a run time error.
   if(length(theta4) ~= length(theta2))
       disp("Matrices' lengths not equal");
       quit(1);
   end
   A = [];
   b = [];
   % rows in both the matrices are added in a loop with columns in each
   for i = 1:length(theta2)
       temp1 = [(cosd(theta4(i))) (-1 * (cosd(theta2(i)))) (1)];
       temp2 = cosd(theta2(i)-theta4(i));
       A = [A; temp1];
       b = [b; temp2];
   end
   link_ratios = lsqr(A,b);% lsqr applies the least square method to solve the matrix for the 3 unknowns
end
function [a,b,c,d] = get_link_lengths(link_ratios)
% get_link_lengths uses the link ratios and the fixed link to find the
% lengths of the other links.
% Lengths can be negatives therefore absolutes of the calculated values are
% sorted.
    d = 180;
    a = abs(d/link_ratios(1));
    c = abs(d/link_ratios(2));
    b = abs(sqrt(a^2  + c^2 + d^2 -(link_ratios(3) * 2 * a * c)));
end
function transmission_angles = get_transmission_angles(a,b,c,d,lower_limit, upper_limit, steps)
% Transmission angles are calculated using the obtained link lengths and
% and the respective input angles 
    transmission_angles = zeros(1, ((upper_limit - lower_limit)/steps));% initializing  array with zeros before use.
    j = 1;% loop counter
    for i = lower_limit:steps:upper_limit
        m = acosd(((b^2 + c^2) - (a^2 + d^2) + (2 * a * d * cosd(i))) / (2 * b * c));
        transmission_angles(j) = m;%structural errors
        j = j + 1;
    end
end
function structuralErrors = get_structural_errors(theta2, theta4, link_ratios)
% Structural error is basically the difference between the left side of the
% freudeinsten's equation and the right side.
    structuralErrors = zeros(1,length(theta4));% theta4 can also be used since they are of the same length
    for i = 1:length(theta2)
        er1 = link_ratios(1)*cosd(theta4(i)) - link_ratios(2)*cosd(theta2(i)) + link_ratios(3) - cosd(theta2(i) - theta4(i));
        structuralErrors(i) = er1;
    end
end