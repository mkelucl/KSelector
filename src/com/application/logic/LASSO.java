package com.application.logic;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.FastMath;
import org.ojalgo.access.Access2D.Factory;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import org.ojalgo.optimisation.Expression.RowColumn;

import com.application.io.TextFileWriter;
import com.application.utils.HelperMethods;

import org.ujmp.core.*;

public class LASSO {
	
	
	
	
	private static final double ETA_DECAY_RATE = 0.8;
	private static final double CLOSE_ENOUGH_THRESHOLD = 0.025; // if its within 5% (plus or minus) of target, then we can stop lasso
	public static double tol = 1e-5;
	public static int MAX_ITERATIONS = 500;
	public static int MAX_RECOMMENDED_PREDICTORS = 20000;
	

	
	
	public static int getVariantsRegion(int variantIndex, int[] regions)
	{
		for (int i = 0; i < regions.length; i++) // go through each region like: [100, 255, 512, etc]
		{
			// if the variant's index is LESS than the region's boundary, then it must be within that region
			if(variantIndex < regions[i]) return i;
		}
		
		System.err.println("variant NOT in any region " + variantIndex );
		return -1;
	}
	
	
	
	// Y is a column vector of n x 1 ( each 
	// X is a (non-square) matrix of n x d.... this is [Individuals] x [SNPs]
	// theta is a d x 1 column vector (this EXCLUDES the intercept)
	// n is the number of observations
	// d is the number of predictors
	// regions: where each region ends, in each region, 1 based. IE 100, 150, 300 etc 
	public static double[][] perform_h2GuidedLASSO_regional(double[] betasWithoutIntercept, double[][] X_variants, double[] Y_phenotypes, double[] snph2Contributions, double[] region_h2s, int[] regions, double hsquared, double shrinkageFactor)
	{
		int i;
		if(snph2Contributions.length > MAX_RECOMMENDED_PREDICTORS)
		{
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			System.out.println("WARNING: YOU HAVE "+ snph2Contributions.length +" PREDICTORS, WHICH IS MORE THAN  THE MAX RECOMMENDED (" + MAX_RECOMMENDED_PREDICTORS + ") - expect java to choke on this ");
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}
		
		// if an extra shrinkage factor was specified, then divide each h2 target by this factor (will have the effect of extra shrinking
		if(shrinkageFactor != -1) hsquared = hsquared * shrinkageFactor;
		else hsquared = hsquared * hsquared;
		
		for (i = 0; i < region_h2s.length; i++)  
		{ 
			if(shrinkageFactor != -1) 
			{
				System.out.println("h2 for region["+(i+1)+"] changed " + (region_h2s[i]) + " -> " + (region_h2s[i] * shrinkageFactor));
				region_h2s[i] = region_h2s[i] * shrinkageFactor; 
			}
			else 
			{
				System.out.println("h2 for region["+(i+1)+"] changed " + (region_h2s[i]) + " -> " + (region_h2s[i] * region_h2s[i]));
				region_h2s[i] = region_h2s[i] * region_h2s[i]; 
			}
			
			
		}
	
		
	
		
		// Doug's formula for obtaining the initial Lambda from the heritability estimate: sqrt (2N/hsq) // N = sample size
		double[] lambdas = new double[region_h2s.length]; // get lambdas for all regions
		for (i = 0; i < region_h2s.length; i++) 
		{	
			//lambdas[i] = FastMath.sqrt( (2* X_variants[0].length )  / region_h2s[i]);
			lambdas[i] = 1;
			System.out.println("lambda for region " + (i+1) +" is: " + lambdas[i]);
		}

		int regionNum;
		double[] regions_h2ExplainedSoFar;
		int[] regions_numNonzeros;
		// double[] regions_meanSNPh2;

		
		double lambda = 0; //FastMath.sqrt( (2* X_variants[0].length )  / hsquared);
		
		
		// 0. Normalize X/Y into Z-score
		System.out.println("starting to standardise Phenotypes/Genotypes");
		
		Y_phenotypes = StatUtils.normalize(Y_phenotypes);
		for (i = 0; i < X_variants.length; i++) 
		{
			// check if there are any monomorphic variants:
			if(StatUtils.variance(X_variants[i]) == 0 ) {  System.err.println("cannot normalise monomorphic variant!, its index: " + i + ". You must exclude these at QC. (use PLINK or similar)"); return null;}
			X_variants[i] = StatUtils.normalize(X_variants[i]);
		}
		System.out.println("finished standardising");
		
		
		
		
		double[] variantsLociIndices = new double[X_variants.length];
		for (i = 0; i < X_variants.length; i++) variantsLociIndices[i] = i;
		
		

		
		// create actual ojAlgo matrices
		final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
		final BasicMatrix X = tmpFactory.columns(X_variants);   System.out.println("X_Matrix is: " + X.toString());
		final BasicMatrix y = tmpFactory.columns(Y_phenotypes);   System.out.println("Y_Matrix is: " + y.toString());
		
		final PhysicalStore.Factory<Double, PrimitiveDenseStore> physicalFactory = PrimitiveDenseStore.FACTORY;
		final PhysicalStore theta = physicalFactory.columns(betasWithoutIntercept);
		

		int j;
		PhysicalStore theta_prev;
		int d = (int) X.countColumns(); 

		// 1. Cache data for aj: 
		final BasicMatrix Xt = X.transpose();  System.out.println("X Transposed Matrix is: " + Xt.toString());  // d x n 
		final BasicMatrix XtX = Xt.multiply(X); // d x m  * m x d =  d x d  ... -> the diagonals are the sum of square of each col of original matrix X
		final BasicMatrix XtX2 = XtX.multiply(2); // this is d x d  // twice as this is a derivative ... the diagonals are aj's for each iteration (twice as this is a derivative)
		System.out.println("XtX2 Transposed Matrix is: " + XtX2.toString());
		
		// 2. cache data for cj
		final BasicMatrix XtY = Xt.multiply(y); // d x n  * n x 1 *  -> d x 1 .. -> a column vector, where is component is the cross product of a column in X with the Y vector
		final BasicMatrix XtY2 = XtY.multiply(2); // crossprod (X,Y) = t(X)%*% Y // twice as this is a derivative
		
		
		// convergence
		int converged = 0;
		int iteration = 0;	
		double aj;
		double cj;
		double sumAbsThetaChange;
		double eta = 0; // how fast we should we be changing lambda
		double h2_threshold = 0; // how close should we ne aiming to hit the target heritability
		double total_h2ExplainedSoFar = 0; // how much h2 the putative causals we've found so far add up to across regions
		double difference;
		double meanSNPh2;

		double old_total_h2ExplainedSoFar = 0;
		
		// initialise learning rate decay
		boolean[] lambda_last_increased = new boolean[d]; // if at last convergence lambda was increased
		double[] eta_decay = new double[d];
		Arrays.fill(eta_decay, 1);
		boolean currentLambdaIncreased = false;
		
		// 3. keep optimising until converged or we've ran out of iterations
		while (converged < 1 && iteration < MAX_ITERATIONS)
		{  
			old_total_h2ExplainedSoFar = getTotalH2Sofar(theta,snph2Contributions/*,regions*/);
			System.out.println("current Iteration is: " +iteration + " // h2: " +  (FastMath.round(old_total_h2ExplainedSoFar * 1000) / 1000.0)   + " / " +  (FastMath.round(hsquared * 1000) / 1000.0)  );
			theta_prev = theta.copy(); // overwrite last iteration's results
			
			
			for (j = 0; j < d; j++) // optimize element j of theta, from the 1st to the Dth
			{
				// aj is twice the sum of squares of Jth column of matrix X: this is equivalent to x[,j]' %*% x[,j] (IE multiply matrix transpose with matrix itself
				aj = (double) XtX2.get(j, j); // grab current aj from the XtX2's diagonals... we have cached these in a matrix, as these are not dependent on theta
				
				
				// cj is defined as =
				// cj = 2 * Xj' * (  (y   -   X   *  beta) +   xj   * beta[j] )
				// where  (y   -   X   *  beta) = the residuals
				// xj   * beta[j] = prediction for the J'th column, 
				// (  (y   -   X   *  beta) +   xj   * beta[j] ) // add the above back onto the residuals, to make them the 'residuals except Jth col'
				// 2 * Xj' = this comes from the fact that this was a partial derivative
				
				// this can be rearranged by multiplying each term by (2 * Xj') as:
				// cj = (Xj' * y  * 2)  -  (Xj' *  X * 2  *  beta) + ( beta[j] * 2 * Xj' * xj )
				
				// which can then make use of the matrices we have cached outside of the loop as:
				// (2Xj' * y  * 2)  =  XtY2[j] // this multiplies the X's j column, with the Y, to give a scalar which is the sum of their crossproduct * 2
				// (Xj' *  X  * 2  *  beta) = XtX2[j,] *  beta // multiply X's J col by the X (* 2x * beta)  // I dont think there is any special 'meaning' here this is just a result of rearrangements
				// ( beta[j] * 2 * Xj' * xj ) =  beta[j] * aj // as  2 * Xj' * xj   was how aj was calculated in the first place

				// finally, this can be written then as:
			//  cj =         XtY2[j]        -                  XtX2[j,]               *     theta                +        theta[j]         *  aj
				cj = (double)XtY2.get(j, 0) - (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) + (double)theta.get(j,0)   *  aj;
				//          1 x 1           -                            (1 x d       *     d x 1)  =   1 x 1   +         1 x 1            * 1 x 1	
				

	
				// get current variant's region and its corresponding lambda
				regionNum = getVariantsRegion(j,regions);
				lambda = lambdas[regionNum];
				
				
				// Then decide how current Theta should change:
				// this is the same as the "softmax" function:
				if (cj < -lambda) theta.set(j, 0, (cj + lambda)/aj);
				else if (cj > lambda) theta.set(j, 0, (cj - lambda)/aj);
				else theta.set(j, 0, 0.0);	
			}
			
			iteration++;
			

			
			// 4. Check for convergence
			// consider it converged if total absolute change is less than threshold	
			//sumAbsThetaChange = 0;
			//for(j = 0; j < theta.countRows();j++)  sumAbsThetaChange += FastMath.abs(theta.doubleValue(j,0) - theta_prev.doubleValue(j,0));
			//if( sumAbsThetaChange < tol) converged = 1;
			
			total_h2ExplainedSoFar = getTotalH2Sofar(theta,snph2Contributions/*,regions*/);
			double diffPerc = FastMath.abs(old_total_h2ExplainedSoFar - total_h2ExplainedSoFar) / old_total_h2ExplainedSoFar;
			if(diffPerc < 0.01)  converged = 1;
			
			
			// 5. if converged, check how much of original heritability is found by summing the Thetas
			if(converged == 1)
			{ System.out.println("converged at iteration: " + iteration); // TODO: maybe check convergence PER region???
				
			 	// must reset this..
				regions_h2ExplainedSoFar = new double[regions.length];
				regions_numNonzeros = new int[regions.length];
				//regions_meanSNPh2 = new double[regions.length];
				
				total_h2ExplainedSoFar = 0;
				// numNonzeros = 0;
				
				
				// find out how much h2 the putative causals we've found explain so far
				for(j = 0; j < theta.countRows();j++)
				{

					// if this predictor was deemed causal (IE nonzero) then we add its h2 contribution
					if( theta.doubleValue(j,0) != 0) 
					{  
						regionNum = getVariantsRegion(j,regions);
						regions_h2ExplainedSoFar[regionNum] += snph2Contributions[j];  
						regions_numNonzeros[regionNum]++; 
						
						total_h2ExplainedSoFar += snph2Contributions[j];  // also keep track of the total
						
					}
				}
				
				// check if we hit target overall h2
				double percDiff =  FastMath.abs(total_h2ExplainedSoFar - hsquared) / hsquared  ; // if 0.98 or 1.02
				if(percDiff > CLOSE_ENOUGH_THRESHOLD) // if difference to target i not close enough
				{	
					// if its further than a set threshold, then remove converged flag,
					for(j = 0; j < regions_h2ExplainedSoFar.length;j++)
					{
						difference = FastMath.abs(region_h2s[j] - regions_h2ExplainedSoFar[j]);
						
						// need to establish the 'close enough' threshold... 
						meanSNPh2 =  regions_h2ExplainedSoFar[j] / regions_numNonzeros[j];
						if(Double.isNaN(meanSNPh2)) meanSNPh2 = 0.0;	
						
					//	if( (j+1) == 1) System.out.println("for region ["+(j+1)+"], putative causals ("+regions_numNonzeros[j]+") explain h2: " + (FastMath.round(regions_h2ExplainedSoFar[j] * 1000) / 1000.0) + " out of real h2: " + (FastMath.round(region_h2s[j] * 1000) / 1000.0) + " / with difference: " + (FastMath.round(difference * 1000) / 1000.0) + " / (meanSNPh2: "+meanSNPh2+") // out of total h2: " + (FastMath.round(total_h2ExplainedSoFar * 1000) / 1000.0) );
						
					//	if( (j+1) == 2) System.out.println("for region ["+(j+1)+"], putative causals ("+regions_numNonzeros[j]+") explain h2: " + (FastMath.round(regions_h2ExplainedSoFar[j] * 1000) / 1000.0) + " out of real h2: " + (FastMath.round(region_h2s[j] * 1000) / 1000.0) + " / with difference: " + (FastMath.round(difference * 1000) / 1000.0) + " / (meanSNPh2: "+meanSNPh2+") // out of total h2: " + (FastMath.round(total_h2ExplainedSoFar * 1000) / 1000.0) );
						
						
						
						
						if(difference >= meanSNPh2) // if we are further away from target than the average SNP's h2 in the region
						{
							converged = 0; // remove converged flag (this will be set to 0, if ANY of the regions fail to have converged...
	
							// determine learning rate: this should be proportionate to: % of how far off we are from h2 target
							eta = lambdas[j] * (difference / region_h2s[j]); // IE if we are 10% off target, then learning rate will be 10% of lambda..
							
							if(eta  > lambdas[j] *0.9 ) eta = lambdas[j] *0.9;
							
							eta = eta * eta_decay[j];
							
							//String changing = (regions_h2ExplainedSoFar[j] < region_h2s[j]) ? "decreasing" : "increasing";
							//if( (j+1) == 4) System.out.println("too far off target, we are " + changing + " region ["+(j+1)+"]'s lambda("+(FastMath.round(lambdas[j] * 1000) / 1000.0)+") by eta: " + (FastMath.round(eta * 1000) / 1000.0));
						//	if( (j+1) == 2) System.out.println("too far off target, we are " + changing + " region ["+(j+1)+"]'s lambda("+(FastMath.round(lambdas[j] * 1000) / 1000.0)+") by eta: " + (FastMath.round(eta * 1000) / 1000.0));
							
							
							// if its LESS than target heritability, then we need to be more 'liberal' and consider more predictors -> decrease Lambda
							if(regions_h2ExplainedSoFar[j] < region_h2s[j]) { lambdas[j] -=eta;  currentLambdaIncreased = false; }
								
							// if its MORE than target heritability, then we need to be more 'strict' and consider less predictors -> increase Lambda
							else { lambdas[j] +=eta; currentLambdaIncreased = true;}
							
							// if its after the first iteration, and we have 'overshot' the target, IE we have changed direction from increasing to decreasing lambda
							if( iteration != 1 &&  currentLambdaIncreased != lambda_last_increased[j] )
							{
								eta_decay[j] = eta_decay[j] * ETA_DECAY_RATE;
								//if( (j+1) == 4) System.out.println("for region ["+(j+1)+"], we have changed eta decay to: " + eta_decay[j]);
							}
							
							lambda_last_increased[j] = currentLambdaIncreased;
							
							
						}
						
						
						
						
						
						
					}
				
				}
				if(converged == 1 ){System.out.println("We've found enough h2 (diff is:"+percDiff+" ), fully converged at: " + iteration );}

			}
			
			
			//System.out.println("at iteration: " +  iteration + " / theta_prev is: " + theta_prev.toString());
			//System.out.println("at iteration: " +  iteration + " / theta_current is: " + theta.toString() + " sumAbsThetaChange: " + sumAbsThetaChange + " isConverged: " + (converged == 1));
  		}

		
        // 5. create output as a (jagged) table: 1st row are thetas, 2nd row the SNP IDs and 3rd row have 2 entries, 1st is the number of iterations it took, and the 2nd is if we have converged
		double[] theta_array = new double[(int) theta.countRows()];
		for(j = 0; j < theta_array.length;j++) theta_array[j] = theta.doubleValue(j,0);
		
		
		// need to create an intermediary table as we cannot sort jagged arrays directly
		// double[][] resultsTable_withoutMeta = new double[2][];
		// dont need this
		// resultsTable_withoutMeta[0] = theta_array;
		// resultsTable_withoutMeta[1] = variantsLociIndices;
		 
		 int pruneCounter =0;
		 List<Double> nonZeroResultThetas = new ArrayList<Double>();
		 List<Double> nonZeroResultIndices = new ArrayList<Double>();
		 for (j = 0; j < theta_array.length; j++)
		 {
			 if(theta_array[j] != 0  ) // if it was found to be unrelated ( IE not 0 )
			 {
				 nonZeroResultThetas.add(theta_array[j]);
				 nonZeroResultIndices.add(variantsLociIndices[j]);
				
				 pruneCounter++;
			 }
		 }
		 
		 
		// get stats:
		regions_h2ExplainedSoFar = new double[regions.length];
		double[] regions_numNonzeros_d = new double[regions.length];
		total_h2ExplainedSoFar = getH2Sofar(theta, regions_h2ExplainedSoFar, regions, regions_numNonzeros_d, snph2Contributions);

		System.out.println("LASSO FINISHED at iteration: " + iteration + ", converged "+(converged == 1)+" // out of total h2: " + (FastMath.round(hsquared * 1000) / 1000.0)  + " it has found " + (FastMath.round(total_h2ExplainedSoFar * 1000) / 1000.0));
		System.out.println("h2 guided LASSO Removed "+ (theta_array.length - pruneCounter)+" unrelated predictors: " );
		 
        double[][] resultsTable = new double[6][];
        resultsTable[0] =  nonZeroResultThetas.stream().mapToDouble(Double::doubleValue).toArray();
        resultsTable[1] =  nonZeroResultIndices.stream().mapToDouble(Double::doubleValue).toArray();
        
        
        // diagnostics
        PearsonsCorrelation corr = new PearsonsCorrelation(); 
        resultsTable[2] = new double[7];
        resultsTable[2][0] = iteration;
        resultsTable[2][1] = converged;
        resultsTable[2][2] = theta_array.length ;
        resultsTable[2][3] = pruneCounter;
        resultsTable[2][4] = total_h2ExplainedSoFar;
        resultsTable[2][5] = hsquared;
        resultsTable[2][6] = corr.correlation(region_h2s, regions_h2ExplainedSoFar);
        resultsTable[3] = region_h2s;
        resultsTable[4] = regions_h2ExplainedSoFar;
        resultsTable[5] = regions_numNonzeros_d;
        
        return resultsTable;
	}
	
	
	public static double[][] perform_h2GuidedLASSO_new(double[] betasWithoutIntercept, double[][] X_variants, double[] Y_phenotypes, double[] snph2Contributions, double hsquared)
	{
		int i;
		if(snph2Contributions.length > MAX_RECOMMENDED_PREDICTORS)
		{
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			System.out.println("WARNING: YOU HAVE "+ snph2Contributions.length +" PREDICTORS, WHICH IS MORE THAN  THE MAX RECOMMENDED (" + MAX_RECOMMENDED_PREDICTORS + ") - expect java to choke on this ");
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}
		
		// double[] regions_meanSNPh2;

		
		double lambda = 1; //FastMath.sqrt( (2* X_variants[0].length )  / hsquared);
		System.out.println("lambda all is: " + lambda);
		
		// 0. Normalize X/Y into Z-score
		System.out.println("starting to standardise Phenotypes/Genotypes");
		
		Y_phenotypes = StatUtils.normalize(Y_phenotypes);
		for (i = 0; i < X_variants.length; i++) 
		{
			// check if there are any monomorphic variants:
			if(StatUtils.variance(X_variants[i]) == 0 ) {  System.err.println("cannot normalise monomorphic variant!, its index: " + i + ". You must exclude these at QC. (use PLINK or similar)"); return null;}
			X_variants[i] = StatUtils.normalize(X_variants[i]);
		}
		System.out.println("finished standardising");
		
		
		
		
		double[] variantsLociIndices = new double[X_variants.length];
		for (i = 0; i < X_variants.length; i++) variantsLociIndices[i] = i;
		
		

		
		// create actual ojAlgo matrices
		final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
		final BasicMatrix X = tmpFactory.columns(X_variants);   System.out.println("X_Matrix is: " + X.toString());
		final BasicMatrix y = tmpFactory.columns(Y_phenotypes);   System.out.println("Y_Matrix is: " + y.toString());
		
		final PhysicalStore.Factory<Double, PrimitiveDenseStore> physicalFactory = PrimitiveDenseStore.FACTORY;
		final PhysicalStore theta = physicalFactory.columns(betasWithoutIntercept);
		

		int j;
		PhysicalStore theta_prev;
		int d = (int) X.countColumns(); 

		// 1. Cache data for aj: 
		final BasicMatrix Xt = X.transpose();  System.out.println("X Transposed Matrix is: " + Xt.toString());  // d x n 
		final BasicMatrix XtX = Xt.multiply(X); // d x m  * m x d =  d x d  ... -> the diagonals are the sum of square of each col of original matrix X
		final BasicMatrix XtX2 = XtX.multiply(2); // this is d x d  // twice as this is a derivative ... the diagonals are aj's for each iteration (twice as this is a derivative)
		System.out.println("XtX2 Transposed Matrix is: " + XtX2.toString());
		
		// 2. cache data for cj
		final BasicMatrix XtY = Xt.multiply(y); // d x n  * n x 1 *  -> d x 1 .. -> a column vector, where is component is the cross product of a column in X with the Y vector
		final BasicMatrix XtY2 = XtY.multiply(2); // crossprod (X,Y) = t(X)%*% Y // twice as this is a derivative
		
		
		// convergence
		int converged = 0;
		int iteration = 0;	
		double aj;
		double cj;
		double sumAbsThetaChange;
		double eta = 0; // how fast we should we be changing lambda
		double h2_threshold = 0; // how close should we ne aiming to hit the target heritability
		double total_h2ExplainedSoFar = 0; // how much h2 the putative causals we've found so far add up to across regions
		double difference;
		double meanSNPh2;

		double old_total_h2ExplainedSoFar = 0;
		
		// initialise learning rate decay
		boolean lambda_last_increased = false; // if at last convergence lambda was increased
		double eta_decay = 1;

		boolean currentLambdaIncreased = false;
		
		// 3. keep optimising until converged or we've ran out of iterations
		while (converged < 1 && iteration < MAX_ITERATIONS)
		{  
			old_total_h2ExplainedSoFar = getTotalH2Sofar(theta,snph2Contributions/*,regions*/);
			System.out.println("current Iteration is: " +iteration + " // h2: " +  (FastMath.round(old_total_h2ExplainedSoFar * 1000) / 1000.0)   + " / " +  (FastMath.round(hsquared * 1000) / 1000.0)  );
			theta_prev = theta.copy(); // overwrite last iteration's results
			
			
			for (j = 0; j < d; j++) // optimize element j of theta, from the 1st to the Dth
			{
				// aj is twice the sum of squares of Jth column of matrix X: this is equivalent to x[,j]' %*% x[,j] (IE multiply matrix transpose with matrix itself
				aj = (double) XtX2.get(j, j); // grab current aj from the XtX2's diagonals... we have cached these in a matrix, as these are not dependent on theta
				
				
				// cj is defined as =
				// cj = 2 * Xj' * (  (y   -   X   *  beta) +   xj   * beta[j] )
				// where  (y   -   X   *  beta) = the residuals
				// xj   * beta[j] = prediction for the J'th column, 
				// (  (y   -   X   *  beta) +   xj   * beta[j] ) // add the above back onto the residuals, to make them the 'residuals except Jth col'
				// 2 * Xj' = this comes from the fact that this was a partial derivative
				
				// this can be rearranged by multiplying each term by (2 * Xj') as:
				// cj = (Xj' * y  * 2)  -  (Xj' *  X * 2  *  beta) + ( beta[j] * 2 * Xj' * xj )
				
				// which can then make use of the matrices we have cached outside of the loop as:
				// (2Xj' * y  * 2)  =  XtY2[j] // this multiplies the X's j column, with the Y, to give a scalar which is the sum of their crossproduct * 2
				// (Xj' *  X  * 2  *  beta) = XtX2[j,] *  beta // multiply X's J col by the X (* 2x * beta)  // I dont think there is any special 'meaning' here this is just a result of rearrangements
				// ( beta[j] * 2 * Xj' * xj ) =  beta[j] * aj // as  2 * Xj' * xj   was how aj was calculated in the first place

				// finally, this can be written then as:
			//  cj =         XtY2[j]        -                  XtX2[j,]               *     theta                +        theta[j]         *  aj
				cj = (double)XtY2.get(j, 0) - (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) + (double)theta.get(j,0)   *  aj;
				//          1 x 1           -                            (1 x d       *     d x 1)  =   1 x 1   +         1 x 1            * 1 x 1	
				
	
				
				// Then decide how current Theta should change:
				// this is the same as the "softmax" function:
				if (cj < -lambda) theta.set(j, 0, (cj + lambda)/aj);
				else if (cj > lambda) theta.set(j, 0, (cj - lambda)/aj);
				else theta.set(j, 0, 0.0);	
			}
			
			iteration++;
			

			
			// 4. Check for convergence
			// consider it converged if total absolute change is less than threshold	
			//sumAbsThetaChange = 0;
			//for(j = 0; j < theta.countRows();j++)  sumAbsThetaChange += FastMath.abs(theta.doubleValue(j,0) - theta_prev.doubleValue(j,0));
			//if( sumAbsThetaChange < tol) converged = 1;
			
			total_h2ExplainedSoFar = getTotalH2Sofar(theta,snph2Contributions/*,regions*/);
			double diffPerc = FastMath.abs(old_total_h2ExplainedSoFar - total_h2ExplainedSoFar) / old_total_h2ExplainedSoFar;
			if(diffPerc < 0.01)  converged = 1;
			
			
			// 5. if converged, check how much of original heritability is found by summing the Thetas
			if(converged == 1)
			{ System.out.println("converged at iteration: " + iteration);
				
				total_h2ExplainedSoFar = 0; // must reset this..
				int numNonzeros = 0;
				
				
				// find out how much h2 the putative causals we've found explain so far
				for(j = 0; j < theta.countRows();j++)
				{
					// if this predictor was deemed causal (IE nonzero) then we add its h2 contribution
					if( theta.doubleValue(j,0) != 0) {  total_h2ExplainedSoFar += snph2Contributions[j];  numNonzeros++; }
				}
				
				// need to establish the 'close enough' threshold...
				meanSNPh2 =  total_h2ExplainedSoFar / numNonzeros;
				if(Double.isNaN(meanSNPh2)) meanSNPh2 = 0.0;
			
			 	
				
				// check if we hit target overall h2
				double percDiff =  FastMath.abs(total_h2ExplainedSoFar - hsquared) / hsquared  ; // if 0.98 or 1.02
				if(percDiff > CLOSE_ENOUGH_THRESHOLD) // if difference to target i not close enough
				{	
					difference = FastMath.abs(hsquared - total_h2ExplainedSoFar);	
					

					converged = 0; // remove converged flag (this will be set to 0

					// determine learning rate: this should be proportionate to: % of how far off we are from h2 target
					eta = lambda * (difference / hsquared); // IE if we are 10% off target, then learning rate will be 10% of lambda..
					
					if(eta  > lambda *0.9 ) eta = lambda *0.9;
					
					eta = eta * eta_decay;
					
					//String changing = (regions_h2ExplainedSoFar[j] < region_h2s[j]) ? "decreasing" : "increasing";
					//if( (j+1) == 4) System.out.println("too far off target, we are " + changing + " region ["+(j+1)+"]'s lambda("+(FastMath.round(lambda * 1000) / 1000.0)+") by eta: " + (FastMath.round(eta * 1000) / 1000.0));
					//	if( (j+1) == 2) System.out.println("too far off target, we are " + changing + " region ["+(j+1)+"]'s lambda("+(FastMath.round(lambda * 1000) / 1000.0)+") by eta: " + (FastMath.round(eta * 1000) / 1000.0));
					
					
					// if its LESS than target heritability, then we need to be more 'liberal' and consider more predictors -> decrease Lambda
					if(total_h2ExplainedSoFar < hsquared) { lambda -=eta;  currentLambdaIncreased = false; }
						
					// if its MORE than target heritability, then we need to be more 'strict' and consider less predictors -> increase Lambda
					else { lambda +=eta; currentLambdaIncreased = true;}
					
					// if its after the first iteration, and we have 'overshot' the target, IE we have changed direction from increasing to decreasing lambda
					if( iteration != 1 &&  currentLambdaIncreased != lambda_last_increased )
					{
						eta_decay = eta_decay * ETA_DECAY_RATE;
						//if( (j+1) == 4) System.out.println("for region ["+(j+1)+"], we have changed eta decay to: " + eta_decay[j]);
					}
					
					lambda_last_increased = currentLambdaIncreased;
					
						
					
				}
				if(converged == 1 ){System.out.println("We've found enough h2 (diff is: "+percDiff+" ), fully converged at: " + iteration );}

			}
			
			
			//System.out.println("at iteration: " +  iteration + " / theta_prev is: " + theta_prev.toString());
			//System.out.println("at iteration: " +  iteration + " / theta_current is: " + theta.toString() + " sumAbsThetaChange: " + sumAbsThetaChange + " isConverged: " + (converged == 1));
  		}

		
        // 5. create output as a (jagged) table: 1st row are thetas, 2nd row the SNP IDs and 3rd row have 2 entries, 1st is the number of iterations it took, and the 2nd is if we have converged
		double[] theta_array = new double[(int) theta.countRows()];
		for(j = 0; j < theta_array.length;j++) theta_array[j] = theta.doubleValue(j,0);
		
		
		// need to create an intermediary table as we cannot sort jagged arrays directly
		// double[][] resultsTable_withoutMeta = new double[2][];
		// dont need this
		// resultsTable_withoutMeta[0] = theta_array;
		// resultsTable_withoutMeta[1] = variantsLociIndices;
		 
		 int pruneCounter =0;
		 List<Double> nonZeroResultThetas = new ArrayList<Double>();
		 List<Double> nonZeroResultIndices = new ArrayList<Double>();
		 for (j = 0; j < theta_array.length; j++)
		 {
			 if(theta_array[j] != 0  ) // if it was found to be unrelated ( IE not 0 )
			 {
				 nonZeroResultThetas.add(theta_array[j]);
				 nonZeroResultIndices.add(variantsLociIndices[j]);
				
				 pruneCounter++;
			 }
		 }
		 
		 
		// get stats:


		total_h2ExplainedSoFar = getTotalH2Sofar(theta,snph2Contributions);

		System.out.println("LASSO FINISHED at iteration: " + iteration + ", converged "+(converged == 1)+" // out of total h2: " + (FastMath.round(hsquared * 1000) / 1000.0)  + " it has found " + (FastMath.round(total_h2ExplainedSoFar * 1000) / 1000.0));
		System.out.println("h2 guided LASSO Removed "+ (theta_array.length - pruneCounter)+" unrelated predictors: " );
		 
        double[][] resultsTable = new double[6][];
        resultsTable[0] =  nonZeroResultThetas.stream().mapToDouble(Double::doubleValue).toArray();
        resultsTable[1] =  nonZeroResultIndices.stream().mapToDouble(Double::doubleValue).toArray();
        
        
        // diagnostics
        resultsTable[2] = new double[7];
        resultsTable[2][0] = iteration;
        resultsTable[2][1] = converged;
        resultsTable[2][2] = theta_array.length ;
        resultsTable[2][3] = pruneCounter;
        resultsTable[2][4] = total_h2ExplainedSoFar;
        resultsTable[2][5] = hsquared;
        
        return resultsTable;
	}
	
	
	/**
	 * Produces an output from a performLASSO() operation
	 * 
	 * @param lassoResults: the result of a LASSO produced by performLASSO()
	 */
	public static void output_LASSO(double[][] lassoResults, int precision, List<String> variantNames, String outputName, String regionStemLoc, int[] regionBoundaries )
	{
		List<String> lines = new ArrayList<String>();
		int i =0;
		String variantID;
		String line;
		// if we had any regions,create their lists...
		ArrayList<List<String>> regionLines = new ArrayList<List<String>>();
		List<String> allSNPs = new ArrayList<String>();
		
		if(regionStemLoc != null && regionBoundaries != null)
		{
			for (i = 0; i < regionBoundaries.length; i++) regionLines.add(new ArrayList<String>());
		}
		
		for (i = 0; i < lassoResults[0].length; i++) // go through all variants
		{
			if(lassoResults[0][i] == 0) continue; // do NOT write zero thetas
			
			// if no real Variant names are supplied, use the default 'SNP' template format
			if(variantNames == null) 
			{
				variantID = "SNP"+((int)lassoResults[1][i] +1);
				line = variantID + "\t" + HelperMethods.roundToDecimals(lassoResults[0][i], precision) ;
				
			}
			
			else   // alternatively, use the variant locus indices to look up the real variant names
			{
				variantID = variantNames.get( (int) lassoResults[1][i] );
				if(variantID.equalsIgnoreCase(".")) System.out.println("WARNING: variant " + variantID + " is not unique, its loci index was: " + (int) lassoResults[1][i]);
				
				line = variantID+ "\t" + HelperMethods.roundToDecimals(lassoResults[0][i], precision);
				
			}
			
			lines.add(line + "\t" + (int) lassoResults[1][i]);
			
			allSNPs.add(variantID);
			
			// if there were any regions specified also write the results out per region
			int regionIndex = 0;
			if(regionStemLoc != null && regionBoundaries != null)
			{
				regionIndex = getVariantsRegion((int)lassoResults[1][i], regionBoundaries);
				
				// get that regions list, and add the variant's ID
				regionLines.get(regionIndex).add(variantID);
			}
		}

		if(allSNPs.size() == 0)
		{
			System.out.println("we didn't find ANY variants!... just use the original SNP list");
			allSNPs = variantNames;
			
			// fill up the lines array with dummy results
			for (i = 0; i < variantNames.size(); i++) lines.add(variantNames.get(i) + "\t" + "-1" + "\t" + i);
	
		}
		// write these separetely, just the snps and their effect sizes
		String outputFileLocation = outputName;
		Path outputFilePath = Paths.get(outputFileLocation);
		TextFileWriter.writeFileOutput(lines, outputFilePath);	
		
		
		// write all SNPs as a single extract list
		outputFilePath = Paths.get(outputFileLocation + ".extract");
		TextFileWriter.writeFileOutput(allSNPs, outputFilePath);	
		
		
		lines = new ArrayList<String>();
		lines.add("CONVERGED:" + "\t" + (lassoResults[2][1] ==1)+ "\t" + "ITERATION:" + "\t" + (int)lassoResults[2][0] );
		lines.add("KEPT:" + "\t" + (int)lassoResults[2][3] + " / " + (int)lassoResults[2][2] + " predictors" + "\t" + "total H2 found:" + "\t" + (FastMath.round(lassoResults[2][4] * 1000) / 1000.0) + " / " +  (FastMath.round(lassoResults[2][5] * 1000) / 1000.0)  );
		
		if(lassoResults[0].length != 0 && regionStemLoc != null && regionBoundaries != null)
		{	
			boolean skippedRegions = false;
			
			lines.add("REGION_NUM" + "\t" + "FOUND_REGION_h2" + "\t" + "TARGET_REGION_h2" + "\t" + "DELTA_(r=" + (FastMath.round(lassoResults[2][6] * 1000) / 1000.0) + ")"+ "\t" + "NUM PREDICTORS" );
		
			for (i = 0; i < lassoResults[3].length; i++)
				lines.add( (i+1) + "\t" + (FastMath.round(lassoResults[4][i] * 10000) / 10000.0) + "\t" + (FastMath.round(lassoResults[3][i] * 10000) / 10000.0) + "\t" + (FastMath.round( (lassoResults[3][i] - lassoResults[4][i]) * 10000) / 10000.0) + "\t" + (int)lassoResults[5][i]);


			// write out region SNPs
			int regionCounter = 0;
			for (i = 0; i < regionBoundaries.length; i++) 
			{
				// first delete old region files - we have to do this, as we may skip certain regions, which means that there would be old region files left
				outputFilePath = Paths.get(regionStemLoc + (i+1));
				try {Files.delete(outputFilePath);} catch (IOException e) {System.err.println("could not delete file: " + outputFilePath.toString()); e.printStackTrace();}
				
				if(regionLines.get(i).size() > 0)
				{
					outputFilePath = Paths.get(regionStemLoc + (regionCounter+1));
					TextFileWriter.writeFileOutput(regionLines.get(i), outputFilePath);
		
					System.out.println("written " + regionLines.get(i).size() + " variants to " + outputFilePath.toString());
					regionCounter++;
				} else { skippedRegions = true; System.out.println("skipped region empty region: " + (i+1) ); }
			}	
			
			if(skippedRegions) // if any regions were skipped due to being empty, we need to rewrite the region_number.txtfile
			{
				System.out.println("as we've skipped regions, we need to rewrite the region_number.txt file, with new number of regions: " + regionCounter);
				List<String> regionNumlines = new ArrayList<String>();
				regionNumlines.add(Integer.toString(regionCounter) );
				outputFilePath = Paths.get(regionStemLoc + "_number.txt");
				TextFileWriter.writeFileOutput(regionNumlines, outputFilePath);
			}
		}
		
		
		outputFilePath = Paths.get(outputFileLocation + ".log");
		TextFileWriter.writeFileOutput(lines, outputFilePath);	
		
		

		
	}
	
	
	

	private static double getH2Sofar(PhysicalStore theta, double[] regions_h2ExplainedSoFar, int[] regions, double[] regions_numNonzeros_d, double[] snph2Contributions)
	{	
		double total_h2ExplainedSoFar = 0;
		int regionNum;
		// find out how much h2 the putative causals we've found explain so far
		for(int j = 0; j < theta.countRows();j++)
		{

			// if this predictor was deemed causal (IE nonzero) then we add its h2 contribution
			if( theta.doubleValue(j,0) != 0) 
			{  
				regionNum = getVariantsRegion(j,regions);
				regions_h2ExplainedSoFar[regionNum] += snph2Contributions[j];  
				regions_numNonzeros_d[regionNum]++; 
				
				total_h2ExplainedSoFar += snph2Contributions[j];  // also keep track of the total
				
			}
		}
		
		return total_h2ExplainedSoFar;
	}
	
	
	private static double getTotalH2Sofar(PhysicalStore theta,double[] snph2Contributions/*, int[] regions*/) 
	{

		double total_h2ExplainedSoFar = 0;
		//int regionNum;
		total_h2ExplainedSoFar = 0;
		for(int j = 0; j < theta.countRows();j++)
		{
			
			// if this predictor was deemed causal (IE nonzero) then we add its h2 contribution
			if( theta.doubleValue(j,0) != 0) 
			{  
				//regionNum = getVariantsRegion(j,regions);
				//regions_h2ExplainedSoFar[regionNum] += snph2Contributions[j];  
				//regions_numNonzeros[regionNum]++; 
			//	System.out.println("adding h2 contribution: " + snph2Contributions[j] );
				
				total_h2ExplainedSoFar += snph2Contributions[j];  // also keep track of the total
			}
		}
		//System.out.println("returning h2 contr so far: " + total_h2ExplainedSoFar );
		return total_h2ExplainedSoFar;
	}



	public static double[][] perform_h2GuidedLASSO(double[] betasWithoutIntercept, double[][] X_variants, double[] Y_phenotypes, double[] snph2Contributions, double hsquared, boolean sort, double[] snph2PosContributions)
	{
		
		// Doug's formula for obtaining the initial Lambda from the heritability estimate: sqrt (2N/hsq) // N = sample size
		//double hsquared = 0.5;
		double lambda = FastMath.sqrt( (2* X_variants[0].length )  / hsquared);
		System.out.println("lambda is: " + lambda);
		
		// 0. Normalize X/Y into Z-score
		System.out.println("starting to standardise Phenotypes/Genotypes");
		int i;
		Y_phenotypes = StatUtils.normalize(Y_phenotypes);
		for (i = 0; i < X_variants.length; i++) 
		{
			// check if there are any monomorphic variants:
			if(StatUtils.variance(X_variants[i]) == 0 ) {  System.err.println("cannot normalise monomorphic variant!, its index: " + i + ". You must exclude these at QC. (use PLINK or similar)"); return null;}
			X_variants[i] = StatUtils.normalize(X_variants[i]);
		}
		System.out.println("finished standardising");
		
		
		if(snph2PosContributions !=  null)
		{ System.out.println("per SNP priors were specified for LASSO, applying these");
			multiplySNPPrior(X_variants,snph2PosContributions);
		} else System.out.println("NO priors specified for SNPs, proceeding as normal");
		
		
		
		double[] variantsLociIndices = new double[X_variants.length];
		for (i = 0; i < X_variants.length; i++) variantsLociIndices[i] = i;
		
		

		
		// create actual ojAlgo matrices
		final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
		final BasicMatrix X = tmpFactory.columns(X_variants);   System.out.println("X_Matrix is: " + X.toString());
		final BasicMatrix y = tmpFactory.columns(Y_phenotypes);   System.out.println("Y_Matrix is: " + y.toString());
		
		final PhysicalStore.Factory<Double, PrimitiveDenseStore> physicalFactory = PrimitiveDenseStore.FACTORY;
		final PhysicalStore theta = physicalFactory.columns(betasWithoutIntercept);
		

		int j;
		PhysicalStore theta_prev;
		int d = (int) X.countColumns(); 

		// 1. Cache data for aj: 
		final BasicMatrix Xt = X.transpose();  System.out.println("X Transposed Matrix is: " + Xt.toString());  // d x n 
		final BasicMatrix XtX = Xt.multiply(X); // d x m  * m x d =  d x d  ... -> the diagonals are the sum of square of each col of original matrix X
		final BasicMatrix XtX2 = XtX.multiply(2); // this is d x d  // twice as this is a derivative ... the diagonals are aj's for each iteration (twice as this is a derivative)
		System.out.println("XtX2 Transposed Matrix is: " + XtX2.toString());
		
		// 2. cache data for cj
		final BasicMatrix XtY = Xt.multiply(y); // d x n  * n x 1 *  -> d x 1 .. -> a column vector, where is component is the cross product of a column in X with the Y vector
		final BasicMatrix XtY2 = XtY.multiply(2); // crossprod (X,Y) = t(X)%*% Y // twice as this is a derivative
		
		
		
		
		
	//	double totalSum = 0;
	//	for(int z = 0; z < Xt.getRowsRange(0, 1).count(); z++)
	//	{
	//		double szar =   (double)Xt.getRowsRange(0, 1).get(z) * (double)X.getColumnsRange(0, 1).get(z)  ;
	//		totalSum = totalSum +  szar;
					
	//	}
	//	System.out.println("XtX.get(0): " + XtX.get(0,0) ); 
	//	System.out.println("totalSum (0,0): " + totalSum );
		
	//	 totalSum = 0;
	//	for(int z = 0; z < Xt.getRowsRange(0, 1).count(); z++)
	//	{
	//		double szar =   (double)Xt.getRowsRange(0, 1).get(z) * (double)X.getColumnsRange(902, 903).get(z)  ;
	//		System.out.println("totalSum: " + totalSum + " with multiplying col/row: " + z + " / ("+szar+")" + " // Xt.getRowsRange(0, 1).get(z): " + Xt.getRowsRange(0, 1).get(z)  + " // X.getColumnsRange(902, 903).get(z): " + X.getColumnsRange(902, 903).get(z));
	//		totalSum = totalSum +  szar;
		
				
	//	}
	//		System.out.println("totalSum (0,902): " + totalSum ); // this is NaN ???
	//		System.out.println("XtX.get(0,902): " + XtX.get(0,902) ); // this is NaN ???
	//		System.out.println("X.get(459, 902): " + X.get(459, 902) );
			
		//	System.out.println("X.get(0, 902): " + X.get(0, 902) );
			
	//	System.out.println("XtX.getRowsRange(0, 1).get(902): " + XtX.getRowsRange(0, 1).get(902) ); // this is NaN ???
		//System.out.println("XtX.getRowsRange(0, 1).get(902): " + XtX.getRowsRange(0, 1).get(902) );
		
		// convergence
		int converged = 0;
		int iteration = 0;	
		double aj;
		double cj;
		double sumAbsThetaChange;
		double eta = 0; // how fast we should we be changing lambda
		double h2_threshold = 0; // how close should we ne aiming to hit the target heritability
		double h2ExplainedSoFar = 0; // how much h2 the putative causals we've found so far add up to
		double difference;
		int numNonzeros;
		double meanSNPh2;

		
		// 3. keep optimising until converged or we've ran out of iterations
		while (converged < 1 && iteration < MAX_ITERATIONS)
		{ System.out.println("current iteration: " +iteration);
			theta_prev = theta.copy(); // overwrite last iteration's results
			
			
			for (j = 0; j < d; j++) // optimize element j of theta, from the 1st to the Dth
			{
				// aj is twice the sum of squares of Jth column of matrix X: this is equivalent to x[,j]' %*% x[,j] (IE multiply matrix transpose with matrix itself
				aj = (double) XtX2.get(j, j); // grab current aj from the XtX2's diagonals... we have cached these in a matrix, as these are not dependent on theta
				
				
				// cj is defined as =
				// cj = 2 * Xj' * (  (y   -   X   *  beta) +   xj   * beta[j] )
				// where  (y   -   X   *  beta) = the residuals
				// xj   * beta[j] = prediction for the J'th column, 
				// (  (y   -   X   *  beta) +   xj   * beta[j] ) // add the above back onto the residuals, to make them the 'residuals except Jth col'
				// 2 * Xj' = this comes from the fact that this was a partial derivative
				
				// this can be rearranged by multiplying each term by (2 * Xj') as:
				// cj = (Xj' * y  * 2)  -  (Xj' *  X * 2  *  beta) + ( beta[j] * 2 * Xj' * xj )
				
				// which can then make use of the matrices we have cached outside of the loop as:
				// (2Xj' * y  * 2)  =  XtY2[j] // this multiplies the X's j column, with the Y, to give a scalar which is the sum of their crossproduct * 2
				// (Xj' *  X  * 2  *  beta) = XtX2[j,] *  beta // multiply X's J col by the X (* 2x * beta)  // I dont think there is any special 'meaning' here this is just a result of rearrangements
				// ( beta[j] * 2 * Xj' * xj ) =  beta[j] * aj // as  2 * Xj' * xj   was how aj was calculated in the first place

				// finally, this can be written then as:
			//  cj =         XtY2[j]        -                  XtX2[j,]               *     theta                +        theta[j]         *  aj
				cj = (double)XtY2.get(j, 0) - (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) + (double)theta.get(j,0)   *  aj;
				//          1 x 1           -                            (1 x d       *     d x 1)  =   1 x 1   +         1 x 1            * 1 x 1	
				

			/*	if( j == 0)
				{
					System.out.println("at iteration " + iteration+ "  // cj is: " + cj + " lambda is: " + lambda + " aj is: " + aj);
					
					//System.out.println("at iteration " + iteration+ "  // (double)XtY2.get(j, 0) is: " + (double)XtY2.get(j, 0) + " (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) is: " + (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) + " (double)theta.get(j,0)   *  aj is: " + (double)theta.get(j,0)   *  aj);
					
					// System.out.println("at iteration " + iteration  + " XtX2.getRowsRange(j, j+1).multiply(theta) is: " + XtX2.getRowsRange(j, j+1).multiply(theta) + " / theta: " + theta + " / XtX2.getRowsRange(j, j+1): " + XtX2.getRowsRange(j, j+1) );
					
					System.out.println("at iteration " + iteration  + " XtX2.getRowsRange(j, j+1).multiply(theta) is: " + XtX2.getRowsRange(j, j+1).multiply(theta) + " / theta.get(0, 0): " + theta.get(0, 0) + " / XtX2.getRowsRange(j, j+1).get(0,0): " + XtX2.getRowsRange(j, j+1).get(0,0) );
					for(int z = 0; z < XtX2.getRowsRange(j, j+1).count(); z++)
					{
						if(Double.isNaN( (double)XtX2.getRowsRange(j, j+1).get(z) )  ) System.out.println("XtX2.getRowsRange(j, j+1).get("+z+"): " + XtX2.getRowsRange(j, j+1).get(z) );
					}

				
				}*/
				
				// Then decide how current Theta should change:
				// this is the same as the "softmax" function:
				if (cj < -lambda) theta.set(j, 0, (cj + lambda)/aj);
				else if (cj > lambda) theta.set(j, 0, (cj - lambda)/aj);
				else theta.set(j, 0, 0.0);	
			}
			
			iteration++;
			
			// 4. Check for convergence
			// consider it converged if total absolute change is less than threshold	
			sumAbsThetaChange = 0;
			for(j = 0; j < theta.countRows();j++)  sumAbsThetaChange += FastMath.abs(theta.doubleValue(j,0) - theta_prev.doubleValue(j,0));
			if( sumAbsThetaChange < tol) converged = 1;
			
			// 5. if converged, check how much of original heritability is found by summing the Thetas
			if(converged == 1)
			{ System.out.println("CONVERGED AT ITERATION: " + iteration);
				h2ExplainedSoFar = 0; // must reset this..
				numNonzeros = 0;
				
				
				// find out how much h2 the putative causals we've found explain so far
				for(j = 0; j < theta.countRows();j++)
				{
					// if this predictor was deemed causal (IE nonzero) then we add its h2 contribution
					if( theta.doubleValue(j,0) != 0) {  h2ExplainedSoFar += snph2Contributions[j];  numNonzeros++; }
				}
				
				// if its further than a set threshold, then remove converged flag,
				difference = FastMath.abs(hsquared - h2ExplainedSoFar);
				
				// need to establish the 'close enough' threshold... for now on
				meanSNPh2 =  h2ExplainedSoFar / numNonzeros;
				if(Double.isNaN(meanSNPh2)) meanSNPh2 = 0.0;
					
				System.out.println("putative causals ("+numNonzeros+") explain h2: " + (FastMath.round(h2ExplainedSoFar * 1000) / 1000.0) + " out of real h2: " + (FastMath.round(hsquared * 1000) / 1000.0) + " / with difference: " + (FastMath.round(difference * 1000) / 1000.0) + " / (meanSNPh2: "+meanSNPh2+")" );
				
				
				if(difference >= meanSNPh2) // if we are further away from target, than the average SNP's h2
				{
					converged = 0; // remove converged flag
				
					// determine learning rate: this should be proportionate to: % of how far off we are from h2 target
					eta = lambda * (difference / hsquared); // 
					
				//	if(eta  > lambda *0.9 ) eta = lambda *0.9;
					
					
					String changing = (h2ExplainedSoFar < hsquared) ? "decreasing" : "increasing";
					System.out.println("too far off target, we are " + changing + " lambda("+(FastMath.round(lambda * 10000) / 10000.0)+") by eta: " + (FastMath.round(eta * 10000) / 10000.0));
					
					
					// if its LESS than target heritability, then we need to be more 'liberal' and consider more predictors -> decrease Lambda
					if(h2ExplainedSoFar < hsquared) lambda -=eta;
						
					// if its MORE than target heritability, then we need to be more 'strict' and consider less predictors -> increase Lambda
					else lambda +=eta;
					
				}
				else {System.out.println("We've found enough h2, fully converged at: " + iteration);}


			}
			
			
			//System.out.println("at iteration: " +  iteration + " / theta_prev is: " + theta_prev.toString());
			//System.out.println("at iteration: " +  iteration + " / theta_current is: " + theta.toString() + " sumAbsThetaChange: " + sumAbsThetaChange + " isConverged: " + (converged == 1));
  		}
	  
		
		
		
		
        // 5. create output as a (jagged) table: 1st row are thetas, 2nd row the SNP IDs and 3rd row have 2 entries, 1st is the number of iterations it took, and the 2nd is if we have converged
		double[] theta_array = new double[(int) theta.countRows()];
		for(j = 0; j < theta_array.length;j++) theta_array[j] = theta.doubleValue(j,0);
		
		
		// need to create an intermediary table as we cannot sort jagged arrays directly
		 double[][] resultsTable_withoutMeta = new double[2][];
		 resultsTable_withoutMeta[0] = theta_array;
		 resultsTable_withoutMeta[1] = variantsLociIndices;
		 
		 // prune negative
		/* int pruneCounter =0;
		 for (j = 0; j < resultsTable_withoutMeta[0].length; j++)
		 {
			 if(snph2Contributions[j] <= 0  && resultsTable_withoutMeta[0][j] != 0 ) { resultsTable_withoutMeta[0][j] = 0;  pruneCounter++;}
			
		 }
		 System.out.println("number of variants pruned: " + pruneCounter);
		*/
		// sort table by the P values so that the most significant results will come first
        if(sort) resultsTable_withoutMeta = HelperMethods.sortTableByRow(resultsTable_withoutMeta,0,true); // sort them on 1st row, as thats where Effect sizes are

		 
		 
        double[][] resultsTable = new double[3][];
        resultsTable[0] =  resultsTable_withoutMeta[0];
        resultsTable[1] =  resultsTable_withoutMeta[1];
        resultsTable[2] = new double[2];
        resultsTable[2][0] = iteration;
        resultsTable[2][1] = converged;

        
        
        
        return resultsTable;
	}

	
	// PLINK's lasso lambda: it calculates both a minimum and a maximum
	// #define NLAMBDA 100 -> magic number for lambda
	// double sqrt_n_recip = sqrt(1.0 / num_samples);
	// sige = sqrt(1.0 - lasso_h2 + 1.0 / num_samples); // sige is sqrt( 1- the heritability, + 1 over the number of sample count )
	// zz = sige * sqrt_n_recip;
	
	// the minimum lambda is then calculated as such:
	// lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS) * zz // the median of some array multiplied by the 'zz' above... WTF
	// 
	
	// Then they calculate a Max lambda as well...
	
	
	// then they go through a 100 iterations of using different lambdas within the calculated max and min
	// and then somehow calculate an error (not sure how), and then I suppose they choose the best one
	// lambda = exp(loghi - logdelta * ((double)((int32_t)lambi)));
	
	
	
	
			//Lasso via coordinate descent: the "Shooting Algorithm" of Fu (1998). Adapted from pseudocode algorithm 13.1 of Murphy (2012) and matlab code LassoShooting.m by Mark Schmidt.
	
			// Y is a column vector of n x 1 ( each 
			// X is a (non-square) matrix of n x d.... this is [Individuals] x [SNPs]
			// theta is a d x 1 column vector (this EXCLUDES the intercept)
			// n is the number of observations
			// d is the number of predictors
			public static double[][] performLASSO(double[] betasWithoutIntercept, double[][] X_variants, double[] Y_phenotypes, double lambda, boolean sort)
			{  		
				
				// Doug's formula for obtaining the initial Lambda from the heritability estimate: sqrt (2N/hsq) // N = sample size
				double hsquared = 0.5;
				//lambda = FastMath.sqrt( (2* X_variants[0].length )  / hsquared);
				System.out.println("lambda is: " + lambda);
				// 0. Normalize X/Y into Z-score
				int i;
				Y_phenotypes = StatUtils.normalize(Y_phenotypes);
				for (i = 0; i < X_variants.length; i++) X_variants[i] = StatUtils.normalize(X_variants[i]);
				
				
				double[] variantsLociIndices = new double[X_variants.length];
				for (i = 0; i < X_variants.length; i++) variantsLociIndices[i] = i;
				
				
				// create actual ojAlgo matrices
				final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
				final BasicMatrix X = tmpFactory.columns(X_variants);  // System.out.println("X_Matrix is: " + X_Matrix.toString());
				final BasicMatrix y = tmpFactory.columns(Y_phenotypes);  // System.out.println("Y_Matrix is: " + Y_Matrix.toString());
				
				final PhysicalStore.Factory<Double, PrimitiveDenseStore> physicalFactory = PrimitiveDenseStore.FACTORY;
				final PhysicalStore theta = physicalFactory.columns(betasWithoutIntercept);
				

				int j;
				PhysicalStore theta_prev;
				int d = (int) X.countColumns(); 

				// 1. Cache data for aj: 
				final BasicMatrix Xt = X.transpose(); // d x n
				final BasicMatrix XtX = Xt.multiply(X); // d x m  * m x d =  d x d  ... -> the diagonals are the sum of square of each col of original matrix X
				final BasicMatrix XtX2 = XtX.multiply(2); // this is d x d  // twice as this is a derivative ... the diagonals are aj's for each iteration (twice as this is a derivative)
				
				// 2. cache data for cj
				final BasicMatrix XtY = Xt.multiply(y); // d x n  * n x 1 *  -> d x 1 .. -> a column vector, where is component is the cross product of a column in X with the Y vector
				final BasicMatrix XtY2 = XtY.multiply(2); // crossprod (X,Y) = t(X)%*% Y // twice as this is a derivative


				int converged = 0;
				int iteration = 0;	
				double aj;
				double cj;
				double sumAbsThetaChange;

				// 3. keep optimising until converged or we've ran out of iterations
				while (converged < 1 && iteration < MAX_ITERATIONS)
				{ System.out.println("current Iteration is: " +iteration);
					theta_prev = theta.copy(); // overwrite last iteration's results
					
					
					for (j = 0; j < d; j++) // optimize element j of theta, from the 1st to the Dth
					{
						// aj is twice the sum of squares of Jth column of matrix X: this is equivalent to x[,j]' %*% x[,j] (IE multiply matrix transpose with matrix itself
						aj = (double) XtX2.get(j, j); // grab current aj from the XtX2's diagonals... we have cached these in a matrix, as these are not dependent on theta
						
						
						// cj is defined as =
						// cj = 2 * Xj' * (  (y   -   X   *  beta) +   xj   * beta[j] )
						// where  (y   -   X   *  beta) = the residuals
						// xj   * beta[j] = prediction for the J'th column, 
						// (  (y   -   X   *  beta) +   xj   * beta[j] ) // add the above back onto the residuals, to make them the 'residuals except Jth col'
						// 2 * Xj' = this comes from the fact that this was a partial derivative
						
						// this can be rearranged by multiplying each term by (2 * Xj') as:
						// cj = (Xj' * y  * 2)  -  (Xj' *  X * 2  *  beta) + ( beta[j] * 2 * Xj' * xj )
						
						// which can then make use of the matrices we have cached outside of the loop as:
						// (2Xj' * y  * 2)  =  XtY2[j] // this multiplies the X's j column, with the Y, to give a scalar which is the sum of their crossproduct * 2
						// (Xj' *  X  * 2  *  beta) = XtX2[j,] *  beta // multiply X's J col by the X (* 2x * beta)  // I dont think there is any special 'meaning' here this is just a result of rearrangements
						// ( beta[j] * 2 * Xj' * xj ) =  beta[j] * aj // as  2 * Xj' * xj   was how aj was calculated in the first place

						// finally, this can be written then as:
					//  cj =         XtY2[j]        -                  XtX2[j,]               *     theta                +        theta[j]         *  aj
						cj = (double)XtY2.get(j, 0) - (double)( (XtX2.getRowsRange(j, j+1).multiply(theta)).get(0, 0) ) + (double)theta.get(j,0)   *  aj;
						//          1 x 1           -                            (1 x d       *     d x 1)  =   1 x 1   +         1 x 1            * 1 x 1	
						

						// Then decide how current Theta should change:
						// this is the same as the "softmax" function:
						if (cj < -lambda) theta.set(j, 0, (cj + lambda)/aj);
						else if (cj > lambda) theta.set(j, 0, (cj - lambda)/aj);
						else theta.set(j, 0, 0.0);	
					}
					
					iteration++;
					
					// 4. Check for convergence
					// consider it converged if total absolute change is less than threshold	
					sumAbsThetaChange = 0;
					for(j = 0; j < theta.countRows();j++)  sumAbsThetaChange += FastMath.abs(theta.doubleValue(j,0) - theta_prev.doubleValue(j,0));
					if( sumAbsThetaChange < tol) converged = 1;
					
					//System.out.println("at iteration: " +  iteration + " / theta_prev is: " + theta_prev.toString());
					//System.out.println("at iteration: " +  iteration + " / theta_current is: " + theta.toString() + " sumAbsThetaChange: " + sumAbsThetaChange + " isConverged: " + (converged == 1));
		  		}
			  
				
				
				
				
		        // 5. create output as a (jagged) table: 1st row are thetas, 2nd row the SNP IDs and 3rd row have 2 entries, 1st is the number of iterations it took, and the 2nd is if we have converged
				double[] theta_array = new double[(int) theta.countRows()];
				for(j = 0; j < theta_array.length;j++) theta_array[j] = theta.doubleValue(j,0);
				
				
				// need to create an intermediary table as we cannot sort jagged arrays directly
				 double[][] resultsTable_withoutMeta = new double[2][];
				 resultsTable_withoutMeta[0] = theta_array;
				 resultsTable_withoutMeta[1] = variantsLociIndices;
				
				// sort table by the P values so that the most significant results will come first
		        if(sort) resultsTable_withoutMeta = HelperMethods.sortTableByRow(resultsTable_withoutMeta,0,true); // sort them on 1st row, as thats where Effect sizes are

				 
				 
		        double[][] resultsTable = new double[3][];
		        resultsTable[0] =  resultsTable_withoutMeta[0];
		        resultsTable[1] =  resultsTable_withoutMeta[1];
		        resultsTable[2] = new double[2];
		        resultsTable[2][0] = iteration;
		        resultsTable[2][1] = converged;

		        
		        
		        
		        return resultsTable;
			}
			
			
			

			
			
			
			public static List<String> load_LASSO_RESULTS(String causalVarsFilename)
			{
				String separatorPattern = "\t";
				
				
				List<String> causalVariants = new ArrayList<String>();
				int counter = 0;
				Scanner s;
				String nextItem;
				String[] nextItemArray;
				
				int colOffset = causalVarsFilename.contains(".lasso") ? 1 : 0; // depending on if results file was generated by KSelector, or PLINK we will use an array offset of either (KSelector)0 or 1 (PLINK)
				int rowOffset = causalVarsFilename.contains(".lasso") ? 0 : 1; // need to skip different number of rows for the header...
				
				String variantString;
				try {
					s = new Scanner (new BufferedReader( new FileReader(causalVarsFilename)));
					//s.useDelimiter(System.getProperty("line.separator")); // OS independent line delimiter
					while(s.hasNext() ) // keep going through the file
					{
						// PLINK format:
						// CHR	SNP	A1	EFFECT
						// 1	SNP107	A	-0.032637
						
						// KSelect format:
						// CONVERGED:	true	ITERATION:	50
						// SNPID	THETA
						// SNP2019	0.1236
						
						nextItem = s.nextLine(); // 1	SNP13	0.435041 (A B 0.415793 0.002479) // for LDAK 5

						if(counter <= rowOffset){ counter++; continue; } // skip the first lines which is the header
						
						nextItemArray = nextItem.split(separatorPattern); 
						
						 causalVariants.add(nextItemArray[colOffset]);
							

						counter++;
					}
				} catch (FileNotFoundException e) {System.out.println(causalVarsFilename); e.printStackTrace(); }
				
				
				return causalVariants;

			}
			
			
			
			
			
			

			public static void multiplySNPPrior(double[][] X_variants, double[] snpPositiveH2Contr)
			{

				int i;
				
			/*	// I. first we want to normalize the 0,1,2 
				for (i = 0; i < X_variants.length; i++) 
				{
					// check if there are any monomorphic variants:
					if(StatUtils.variance(X_variants[i]) == 0 ) {  System.err.println("cannot normalise monomorphic variant!, its index: " + i + ". You must exclude these at QC. (use PLINK or similar)"); return ;}
					X_variants[i] = StatUtils.normalize(X_variants[i]);
				}
			
				double [] X_variants_1_0_Z = new double[X_variants[1].length];
				for (int j = 0; j < X_variants[1].length; j++)	X_variants_1_0_Z[j] = X_variants[1][j];
				*/	
					
				// calculate the SNP priors: these are the % deviation from the mean
				double[] snpPriors = new double[snpPositiveH2Contr.length];
				double meanSNPh2 = StatUtils.mean(snpPositiveH2Contr); // get mean expected h2 contribution of each SNP
				double SNPi_Prior;
				for (i = 0; i < X_variants.length; i++) // go through each SNP
				{	
					SNPi_Prior = snpPositiveH2Contr[i] / meanSNPh2; // get weight: current / expected
					snpPriors[i] = SNPi_Prior;
					

					// multiply each Individuals genotype value by this prior
					for (int j = 0; j < X_variants[i].length; j++) // go through each individual, and multiply their genotype values by the prior for SNP i
					{
						//if(i == 0 && j == 0) System.out.println("multiplying " + X_variants[i][j] + " by " + SNPi_Prior);
						X_variants[i][j] = X_variants[i][j] * SNPi_Prior;
						//if(i == 0 && j == 0) System.out.println("X_variants[i][j]: " + X_variants[i][j]);
					}
				}
				
				//// DEBUG
				//System.out.println("X_variants[0][0] After Weighting is: " + X_variants[0][0] + " its SNP prior was: " + snpPriors[0]);		
				//System.out.println("X_variants[1][0] After Weighting is: " + X_variants[1][0] + " its SNP prior was: " + snpPriors[1]);
				// DEBUG

				
				// do NOT restandardise it again... as that would reset the priors ( as multiplying a vector by a constant will not alter their Z-scores)

				//double [] X_variants_1_0_Z_AGAIN = X_variants[1];
				
				//System.out.println("finished applying priors, meanSNPh2 was: " + meanSNPh2);
				
				System.out.println("finished applying priors, meanSNPh2 was: " + meanSNPh2);
			}
			
			
}
