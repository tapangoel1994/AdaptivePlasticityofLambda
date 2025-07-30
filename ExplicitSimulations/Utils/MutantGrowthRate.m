% Function reads in the timevector and density vector of resident and
% mutant populations (assuming there is only a resident and a mutant with
% max MOI = 2). It finds peaks of mutant fraction and fits an exponential
% to it to obtain the growth rate.

function [rate,prefactor] = MutantGrowthRate(T,Y)

    [~,c] = size(Y);
    if c > 11
        error('Data is incompatible with analysis assumptions');
    else

        % obtain densities of resident and mutant genomes
        resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10);
        mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11);
        
        % obtain mutant fraction
        mutantfraction = mutant./(resident+mutant);
        
        % find peaks of mutant fraction
        [pks, locs] = findpeaks(mutantfraction,T);
        
        % fit exponential to peaks of mutant fraction
        f1 = fit(locs,pks,'exp1');
        
        
        % return exponential decay rate
        rate = f1.b;
        prefactor = f1.a;
    end
end