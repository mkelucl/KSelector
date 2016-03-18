package com.application.logic;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;
import org.ojalgo.access.Access2D.Factory;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;



import com.application.io.TextFileWriter;
import com.application.utils.HelperMethods;

public class OLS {
	
	
	
	/**
	 * Performs Multiple regression on each variant and returning their associated effect sizes/P-values
	 * NOTE: cant use this for sparse models where n < p, as there will be negative degrees of freedom and model will blow up
	 * @param predictor_X: the predictors ( Genotype matrix, with following signature [variant][samples] )
	 * @param outcome_Y: the outcome (Y)
	 * @param bonferroniCorrect: if the P values should be bonferroni corrected
	 * @param storeAsMinLogP: if P-values should be stored as minLog10 
	 * @param sort: if the results table should be sorted so that most significant results come first
	 * @return a 2d matrix with a signature of [0] = P-values, [1] = Betas, [2] = variants's LociIndices, [3] = variants's SE
	 */
		public static double[][] performMultiRegOLS(double[][] predictor_X, double[] outcome_Y, boolean bonferroniCorrect,  boolean storeAsMinLogP, boolean sort) 
		{
			int i = 0;
			
			
			double[] variantsLociIndices = new double[predictor_X.length];
			for (i = 0; i < predictor_X.length; i++) variantsLociIndices[i] = i;

			// 2. produce Data in standard format: X: 1st col being all 1s with all the non predictor variables , and Y being the outcome vars
			// create all 1 column for intercept:
			double[] interceptBuffer = new double[predictor_X[0].length];
			Arrays.fill(interceptBuffer, 1);
		
			// create data for X input matrix
			double[][] X_2DArray = new double[predictor_X.length + 1][]; // +1 for the all 1s column for the intercept
			X_2DArray[0] = interceptBuffer;
			for(i = 0; i < predictor_X.length; i++) X_2DArray[i+1] = predictor_X[i]; // 
		
			// create data for Y output matrix
			double[][] Y_2DArray = new double[1][];
			Y_2DArray[0] = outcome_Y;
			
			
			// 3. create actual ojAlgo matrices
			final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
			final BasicMatrix X_Matrix = tmpFactory.columns(X_2DArray);  // System.out.println("X_Matrix is: " + X_Matrix.toString());
			final BasicMatrix Y_Matrix = tmpFactory.columns(Y_2DArray);  // System.out.println("Y_Matrix is: " + Y_Matrix.toString());
			

			// 4. get MLE Theta estimates from formula: ThetaML = ( Xt * X )^-1 * Xt * Y
			// ( Xt * X )^-1
			final BasicMatrix Xtransposed_Matrix = X_Matrix.transpose();
			final BasicMatrix XtransposedX_Matrix = Xtransposed_Matrix.multiply(X_Matrix);
			final BasicMatrix XtransposedXInverted_Matrix = XtransposedX_Matrix.invert();
			
			final BasicMatrix XtransposedXY_Matrix = Xtransposed_Matrix.multiply(Y_Matrix); //  Xt * Y
			final BasicMatrix ThetaMLE = XtransposedXInverted_Matrix.multiply(XtransposedXY_Matrix); // ( Xt * X )^-1 * Xt * Y
			
			System.out.println("ThetaMLE is: " + ThetaMLE.toString());
			List<Double> ThetasAsList = (List<Double>) ThetaMLE.toListOfElements();
			
			
			// _______________________________

			// 5. get MLE regression variance estimate from formula:  http://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression
			// G^2 = 1/n * (Y - X*Theta)transposed * (Y - X*Theta)
			long degreesOfFreedom = X_Matrix.countRows() -  X_Matrix.countColumns();  System.out.println("degreesOfFreedom is: " + degreesOfFreedom); // the degrees of freedom is the number of Samples minus the number of Predictors 
			double OnePerN = 1.0 / degreesOfFreedom; // 1/n    //System.out.println("OnePerN is: " + OnePerN + " / X_Matrix.countRows(): " + X_Matrix.countRows());
			                       // we are penalizing it as it is an estimate, otherwise it would be just 1.0 / X_Matrix.countRows()
			final BasicMatrix XtimesTheta_matrix = X_Matrix.multiply(ThetaMLE); // X*Theta      //System.out.println("XtimesTheta_matrix is: " + XtimesTheta_matrix );
			final BasicMatrix YminusXTheta_matrix = Y_Matrix.subtract(XtimesTheta_matrix); // (Y - X*Theta) 
			final BasicMatrix YminusXThetaTransposed_matrix = YminusXTheta_matrix.transpose(); // (Y - X*Theta)transposed
			final BasicMatrix residuals_matrix = YminusXThetaTransposed_matrix.multiply(YminusXTheta_matrix); // (Y - X*Theta)transposed * (Y - X*Theta)
			final BasicMatrix covariance_matrix = residuals_matrix.multiply(OnePerN);
			//List<Double> varianceList = (List<Double>) covariance_matrix.toListOfElements();
			double sigmaSquared = ( (List<Double>) (covariance_matrix.toListOfElements()) ).get(0); System.out.println("variance as scalar is: " + sigmaSquared);
			
			
			// Get misc regression parameters: SE, T-value and P-value
			// _______________________________
			// _______________________________

			// Get Variance - Covariance matrix:
			final BasicMatrix varianceCovariance_matrix = XtransposedXInverted_Matrix.multiply(sigmaSquared);  // System.out.println("varianceCovariance_matrix is: " + varianceCovariance_matrix.toString());
			

			// get diagonal elements of varianceCovariance_matrix (no built in optimised function)
			double [] diagonal = new double [(int) varianceCovariance_matrix.countColumns()];
			for(i = 0; i < diagonal.length; i++) diagonal[i] = varianceCovariance_matrix.doubleValue(i, i);

			
			// the SEs of the slope are the square roots of the diagonal of the variance Covariance matrix
			double [] coefSEs = new double[diagonal.length];
			for(i = 0; i < diagonal.length; i++) coefSEs[i] = FastMath.sqrt(diagonal[i]);
			System.out.println("SEs of coefs are: " + Arrays.toString(coefSEs));

			// calculate T and P values
	        double[] tValues = new double[coefSEs.length];
	        for (i = 0; i < coefSEs.length; i++)  tValues[i] = ThetasAsList.get(i) / coefSEs[i];
	        System.out.println("tValues are: " + Arrays.toString(tValues));
	        
	        // get P-values
	        TDistribution tDistribution = new TDistribution(degreesOfFreedom);
	        double[] pValues = new double[coefSEs.length];
	        for (i = 0; i < tValues.length; i++) {
	        	pValues[i] =   2*(  1- tDistribution.cumulativeProbability(FastMath.abs(tValues[i]))); //System.out.println("UNCORRECTED " + pValues[i] + " / for X" + (i -1));
	        	
	        	// apply bonferroni correction only if requested
	        	if(bonferroniCorrect) pValues[i] = FastMath.min(pValues[i] * coefSEs.length, 1); // Bonferroni correct them ( IE multiply the significance level by the number of tests, IE X values) and cap it to 1

	        	if(storeAsMinLogP) pValues[i] = -FastMath.log10(pValues[i]); // store them as minus log10(p)
	           //System.out.println("-Log(P) " + pValues[i] + " / for X" + (i -1));
	        }
	        System.out.println("pValues are: " + Arrays.toString(pValues));
	        
	        
	        // create output as a table: 1st row are coefs, 2nd row are P values, 3rd row SE's
	        double[][] resultsTable = new double[4][];
	        resultsTable[0] = pValues;
	        resultsTable[1] = new double[ThetasAsList.size()];
	        for ( i = 0; i < ThetasAsList.size(); i++) resultsTable[1][i] = ThetasAsList.get(i);
	        resultsTable[2] = variantsLociIndices; // must store this, in order to be able to match back
	        resultsTable[3] = coefSEs;
	        
	        
	        // sort table by the P values so that the most significant results will come first
	       if(sort) resultsTable = HelperMethods.sortTableByRow(resultsTable,0,storeAsMinLogP); // sort them on 2nd row, as thats where P values are, and store it descending IF we are storing them as minlog10, as then the largest will be most significant

	        
	        return resultsTable;
	       
		}
		
		
		

		/**
		 * Performs basic GWAS analysis on a set of Samples/Variants, by performing UNIVARIATE regression on each variant and returning their associated effect sizes/P-values
		 * This is calibrated so that it gives identical results to PLINK2
		 * @param predictor_X: the predictors ( Genotype matrix, with following signature [variant][samples] )
		 * @param outcome_Y: the outcome (Y)
		 * @param bonferroniCorrect: if the P values should be bonferroni corrected
		 * @param storeAsMinLogP: if P-values should be stored as minLog10 
		 * @param sort: if the results table should be sorted so that most significant results come first
		 * @return a 2d matrix with a signature of [0] = P-values, [1] = Betas, [2] = variants's LociIndices, [3] = variants's SE
		 */
		public static double[][] performGWAS(double[][] predictor_X, double[] outcome_Y, boolean bonferroniCorrect,  boolean storeAsMinLogP, boolean sort)
		{
			int i;
			int j;
	
			 SimpleRegression lm = new SimpleRegression(); 
			 
			 double[] beta = new double[predictor_X.length];
			 double[] pValues = new double[beta.length];
			 double[] coefSEs = new double[beta.length];
			 double[] variantsLociIndices = new double[beta.length];
			 for (i = 0; i < predictor_X.length; i++) variantsLociIndices[i] = i;
			 
			 double[][] tempData;
			 for (i = 0; i < predictor_X.length; i++) // go through all loci
			 {	
				 // ordinary GWAS, we create a new univeriate regression model for each variant
				lm.clear();
				//lm = new SimpleRegression(); 
				// add data into the Linear model, create 2D array, where 1st row are all the X, and 2nd row are all the Y values
				tempData = new double[2][outcome_Y.length];
				tempData[1] = predictor_X[i];
				tempData[0] = outcome_Y;
				// this is rotated the wrong way, so we need to transpose it
				tempData = HelperMethods.transposeMatrix(tempData); // put SNPs on the rows
				lm.addData(tempData);
			
				
				
				// it is still better than doing this, as this would load data 1 by 1, and recalculate stuff for EACH new datapoint...
			//	for(j =0; j< predictor_X[j].length; j++) 
			//		lm.addData(predictor_X[i][j],outcome_Y[j]);  //create LM model, and load data 1 by 1...
				
				
				
				// get Slope (ie Beta)
				beta[i]     = lm.getSlope();	
				pValues[i]  = lm.getSignificance(); // get P value
				coefSEs[i]  = lm.getMeanSquareError(); // store Error

				
	        	// apply bonferroni correction only if requested
	        	if(bonferroniCorrect) pValues[i] = FastMath.min(pValues[i] * beta.length, 1); // Bonferroni correct them ( IE multiply the significance level by the number of tests, IE X values) and cap it to 1
	        															// as Java's 'x.length' refers to the NUMBER OF POINTS, and the number of tests is in the x[0].length ... +1 the intercept... IE 76 Variants, with 100 datapoints, = 77 tests
	        	if(storeAsMinLogP) pValues[i] = -FastMath.log10(pValues[i]); // store them as minus log10P
	          // System.out.println("for SNP("+allVariants.get(variantIndex).START+"): " + variantIndex + " -> -Log(P): " + pValues[variantIndex] + " / P-value: "+lm.getSignificance() );
			 }	

			
			// IV. return results of genetic effects (Betas) and their Significance ( associated P-values)
	        double[][] matrixResult = new double[4][];
	        matrixResult[0] = pValues; 
	        matrixResult[1] = beta;
	        matrixResult[2] = variantsLociIndices; // must store this, in order to be able to match back
	        matrixResult[3] = coefSEs; 

	        
	        
	     // sort table by the P values so that the most significant results will come first
	       if(sort) matrixResult = HelperMethods.sortTableByRow(matrixResult,0,storeAsMinLogP); // sort them on 2nd row, as thats where P values are, and store it descending IF we are storing them as minlog10, as then the largest will be most significant

	        
			return matrixResult;
		}
		
		
		/**
		 * Produces an output from a performGWAS() operation
		 * 
		 * @param matrixResult: the result of a GWAS produced by performGWAS()
		 * @param allVariants: all the variants (must have the same order/index as the matrix results)
		 * @return: lines ready to be written into a file
		 */
		public static void output_GWAS(double[][] gwasResults, int precision )
		{
			int j;

			List<String> lines = new ArrayList<String>();
			int i =0;
			String variantID;
			lines.add("SNPID" + "\t" + "P-value" + "\t" + "BETA" );
			for (i = 0; i < gwasResults[2].length; i++) // go through all variants
			{

				lines.add("SNP"+((int)gwasResults[2][i] +1) + "\t" + HelperMethods.roundToDecimals(gwasResults[0][i], precision) + "\t" + HelperMethods.roundToDecimals(gwasResults[1][i],precision) );
				// 85	variant85	c	86	0.1	0.05
			}

			String outputFileLocation = "KGWASresults.txt";
			Path outputFilePath = Paths.get(outputFileLocation);
			TextFileWriter.writeFileOutput(lines, outputFilePath);	
		}

		
}
