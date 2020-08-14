function Transprob = markov_trans(Mu)

	% markov transtion probabilty matrix
	l = length(Mu);
	Transprob = zeros(l);

	for i=1:l
		for j=1:l
			if i==j
				Transprob(i,j) = 0.7;
			else
				Transprob(i,j) = 0.3/(l-1);
			end
		end
	end

end