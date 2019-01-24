This is the progress report as of 1/23/19 for the henric branch of continuous-variable-evolution.

The objective was to program an agent based approximation to the Henrich 2004 model to see if
it could be an example of the nearly neutral model.

Brief description:  Based on the continuous variable model with 1 subpopulation and 1 attribute.

The atribute value of each individual was mutated multiplicatively on every generation.  The mutation 
distribution is normally distributed with std deviation mutstddev and mean mutation_bias (where
mutstddev and mutation_bias are parameters).

The fitness of an individual was the attribute value.  

Selection was either elitist (winner-take-all) to correspond to the Henrich model, or proportional
selection applied to the fitenss to the fit_power power.  As fit_power goes to infinity, this 
selection goes to elitist.

Thus, selection increases the attribute values (which correspons to skill level in the Henrich model),
and mutation tends to decrease the attribute value.

The non-elitist selection is needed to get population variation (coeficient of variation) results.

There is also a renormalize parameter which resets the attribute value to 1 on every generation.

I did not modify exixting contvar files.  Instead, the modified files all start with h or H.


Desired results;
First, fitness should increase with population size because this is what the Henrich model predicts.

Second, results should correspond to the nearly neutral results for Figures 8, 10, 11, 12 of the 
nearly neutral paper.

Results:

Setting fit_power = 12.0  seems to give both results fairly comparable to elitist slection and
reasonable results for Coef Var.

Results are in data/1_18_19 (preliminary), 
1_20_19  (with mutation_bias = 0.95)
1_21_19  (with mutation bias = 0.98)
The Windows xlsx files include 
pivot charts.

The N_mut attribute coef var results  correspond very well to Figure 10.  N6400_25_N_mut_bias095_ngens3_fp12.xlsx

N6400_25_N_mut_bias095_ngens3_fp12.xlsx
The mutstddev fitness mean results show an increase with popsize. Very approximately matches Figure 12.
But the fitness model is different so we shouldn't expect an exact match.

The attribute coef var results approximately match Figure 8.

However, the nearly neutral fit diff results show that most fit differences are neg neutral which justifies
the nearly neutral assumption.  The Henrich results are quite different with most fit diffs being non-neutral
negative.  This perhaps suggests that the nearly neutral model does not apply.



