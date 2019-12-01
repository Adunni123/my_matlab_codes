function MyCG(A, b)

m=length(A); x=MyZeros(m,1); % set x=0 as IC

r = b - A*x;   % Compute first residual
p = r;   
inf_norm=MyZeros(m,1); two_norm=MyZeros(m,1);

for i=1:m
        
	alpha = (trans(p)*r)/(trans(p)*A*p);
    
	x = x + alpha*p;  % update x
    
	r = b - A*x; % update r
    
    inf_norm(i) = max(MyAbs(r));  % obtain infinity norm of residual
    
    two_norm(i) = sqrt( sum(MyAbs(r).^2) );  % obtain 2-norm of residual
        	
    beta = -(trans(p)*A*r)/(trans(p)*A*p);
    
	p = r + beta * p; % update p
end

disp('The CG solution @(0.06,0.04) = ');
% disp(x);
disp(x(12));

figure
plot(inf_norm, 'r');    % plot infinity-norm
xlim([1,19]); xlabel('Iterations','FontSize',12,'FontWeight','bold'); 
ylabel('\infty - norm', 'FontSize',12,'FontWeight','bold')

figure
plot(two_norm, 'b');    % plot 2-norm
xlim([1,19]); xlabel('Iterations', 'FontSize',12,'FontWeight','bold'); 
ylabel('Two - norm', 'FontSize',12,'FontWeight','bold')

end
