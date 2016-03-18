package com.application.logic;

import java.util.ArrayList;
import java.util.List;

public class PerformanceSummary {


	
	/**
	 * Compares the putative list of variants against a real list
	 * @param putativeCausalVariants: the list of variants that the analysis are believed to be causal
	 * @param howManyToCheck: only this many of the putative causal variants are checked against ( must always be < putativeCausalVariants.size()
	 * @param solutionsCausalVariants: the actual list of real Causals' locus indices
	 * @return: ready to be written contents for a text file
	 */
	public static List<String> compareAgainstSolutions(List<String> putativeCausalVariants, List<String> solutionsCausalVariants, List<String> variantNames, double[] snph2contr, double h2)
	{
		// go through list of solutions
		String variant;
		List<String> lines = new ArrayList<String>();
		int realCausalFoundCounter = 0;
		lines.add("SNPID" + "\t" + "FOUND"  );
		StringBuilder sb;
		int i;
		int positionOfRealCausalOnPutativeList = 0;
		//  
		// I. Go through the real list of causal variants: determine how many were found
		for(i = 0; i < solutionsCausalVariants.size(); i++)
		{
			sb = new StringBuilder();
			
			variant = solutionsCausalVariants.get(i);
			sb.append( variant + "\t");
			
			positionOfRealCausalOnPutativeList = putativeCausalVariants.indexOf(variant);
			if(positionOfRealCausalOnPutativeList != -1) // IE, if its ON the list	
			{
				realCausalFoundCounter++;
				sb.append("YES");
			}
			else sb.append("NO");
			
			lines.add(sb.toString());
		}
		//lines.add("_______________________________________" );
		
		double totalH2Explained = 0.0;
		// II. Go through the putative found variants, and see if there are any that were NOT on the real list, to determine any false positives
		int falsePositiveCounter = 0;
		for(i = 0; i < putativeCausalVariants.size(); i++) // only check the Top hits
		{
			variant = putativeCausalVariants.get(i);
			// if it was NOT on the real list,   then we have got a false positive
			if(solutionsCausalVariants.indexOf(variant) == -1) falsePositiveCounter++;
			
			// also find out how much heritability the putative causals have explained
			totalH2Explained += snph2contr[variantNames.indexOf(variant)];
				
		}
		double FCV =realCausalFoundCounter/(double)solutionsCausalVariants.size();
		double FPR =falsePositiveCounter/(double)putativeCausalVariants.size();
		String finalAnalysis = "Found causal variants: " + realCausalFoundCounter + " / " + solutionsCausalVariants.size() + " ("+Math.round((FCV) * 100.0)+"%)";
		String falsePositives = "False Positives: " + falsePositiveCounter + " / " + putativeCausalVariants.size() + " ("+Math.round((FPR) * 100.0)+"%)";
		
		System.out.println(finalAnalysis);
		System.out.println(falsePositives);
		lines.add(0,finalAnalysis);
		lines.add(0,falsePositives);
		lines.add(0,Math.round((FCV) * 100)/100.0 + "\t" + Math.round((FPR) * 100)/100.0+ "\t" + Math.round((totalH2Explained) * 100)/100.0+"/"+ (Math.round(h2 * 100)/100.0) );
		lines.add(0,"Found_Causal_Variants" + "\t" + "False_Positive_Rate"+ "\t" + "total_h2_explained" );
		
		
		return lines;
	}
}
