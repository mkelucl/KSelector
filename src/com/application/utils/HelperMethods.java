// Copyright (c) 2015, Marton Kelemen
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


package com.application.utils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;



//import com.application.data.StaticData;

public class HelperMethods {
	/**
	 * this is used as a 'safe' method for parsing GT fields as normally these are like this '0|1', however when a call could not be made
	 *  it has a '.' instead which breaks the Java parseInt method...
	 *  so in those cases the call is going to be mapped to '-1' 
	 */
	
	//public static Boolean FORCE_UNIX_FORMAT = false;
	
	public static int tryParseInt(String value) { 

		 try {  
		     return Integer.parseInt( value.trim());  
		  } catch(NumberFormatException nfe) {  
			 
			  
	  
		// if the above failed the first time around, then try to clean the string from non numeric characters
		// as some people think its OK to enter CHROM 16 as 'chrom16'
		// source: http://stackoverflow.com/questions/13440506/java-find-all-numbers-in-the-string-need-check    
		value = value.replaceAll("[^0-9]+", "");
		value = value.trim();
		
			  
		 try {  
		     return Integer.parseInt( value.trim());  
		  } catch(NumberFormatException nfe2) {  
  
					  
		      return -1; //StaticData.NOT_SET;
		  } 
		
		  }
		}
	
	public static Double tryParseDouble(String value) {

		try{
		   return Double.parseDouble(value.trim());
		}
		catch(NumberFormatException nfe){
		   // System.out.println("The string is not formatted correctly");
			return Double.NaN; // StaticData.NOT_SET_D;
		}
	}


	
	public static String[] concat(String[] a, String[] b) {
		   int aLen = a.length;
		   int bLen = b.length;
		   String[] c= new String[aLen+bLen];
		   System.arraycopy(a, 0, c, 0, aLen);
		   System.arraycopy(b, 0, c, aLen, bLen);
		   return c;
		}
	
	
	public static Double roundToDecimals(Double number, int places)
	{
		if(places == -1) return number; // if there is no rounding to take place
		
		Double roundingFactor = Math.pow(10, places);
		return Math.round(number * roundingFactor) / roundingFactor;
		
	}
	
	// http://stackoverflow.com/questions/1564832/how-do-i-do-a-deep-copy-of-a-2d-array-in-java
	public static byte[][] deepCopy(byte[][] original) {
	    if (original == null) {
	        return null;
	    }

	    final byte[][] result = new byte[original.length][];
	    for (int i = 0; i < original.length; i++) {
	        result[i] = Arrays.copyOf(original[i], original[i].length);
	        // For Java versions prior to Java 6 use the next:
	        // System.arraycopy(original[i], 0, result[i], 0, original[i].length);
	    }
	    return result;
	}
	
	
	
	
	
	
	// checks if an Array contains an element: most efficient for small arrays: http://www.programcreek.com/2014/04/check-if-array-contains-a-value-java/
	public static boolean arrayContains(byte[][] arr, byte[] targetValue) {
		for(byte[] s: arr){
			if(s.equals(targetValue))
				return true;
		}
		return false;
	}
	
	   
		public static ArrayList<byte[]> deepCopyArrayList(ArrayList<byte[]> original) {
		    if (original == null) {
		        return null;
		    }

		    final ArrayList<byte[]> result = new  ArrayList<byte[]>();
		    for (int i = 0; i < original.size(); i++) {
		        result.add(Arrays.copyOf(original.get(i), original.get(i).length));
		        // For Java versions prior to Java 6 use the next:
		        // System.arraycopy(original[i], 0, result[i], 0, original[i].length);
		    }
		    return result;
		}
		/**
		 * Returns a pseudo-random number between min and max, inclusive.
		 * The difference between min and max can be at most
		 * <code>Integer.MAX_VALUE - 1</code>.
		 *
		 * @param min Minimum value
		 * @param max Maximum value.  Must be greater than min.
		 * @return Integer between min and max, inclusive.
		 * @see java.util.Random#nextInt(int)
		 * 
		 *  source: http://stackoverflow.com/questions/363681/generating-random-integers-in-a-range-with-java
		 */
		public static int randInt(int min, int max) {

		    // NOTE: Usually this should be a field rather than a method
		    // variable so that it is not re-seeded every call.
		    Random rand = new Random();
		    
			return randInt(min,max,rand);
		}	


		public static int randInt(int min, int max, Random rand) {
		    // nextInt is normally exclusive of the top value,
		    // so add 1 to make it inclusive
		    int randomNum = rand.nextInt((max - min) + 1) + min;

		    return randomNum;
		}
	// _____________________________________________
		
		// basic stats functions
		
		public static double getSum(double[] X)
		{
			if(X == null || X.length == 0) return 0.0;
			Double total = 0.0;
			for(int i =0; i < X.length; i++) 
				{	
					total+= X[i];
					//if(Double.isNaN(total)) System.out.println("total is NaN at i: " + i);
					if(Double.isNaN(X[i])) System.out.println("X[i] is NaN at i: " + i);
				}
			
			return total;
		}
		
		

		public static Double getMean(double[] X)
		{	
			return getSum(X) / X.length;
		}
		public static Double getSumsOfSquares(double[] X, Double mean)
		{
			Double total = 0.0;
			for(int i =0; i < X.length; i++) total+= ( (X[i] - mean) * (X[i] - mean) );
			
			return total;
		}
		
		public static Double getProductsAboutMean(double[] X, Double meanX,double[] Y, Double meanY)
		{
			Double total = 0.0;
			for(int i =0; i < X.length; i++) total+= ( (X[i] - meanX) * (Y[i] - meanY) );
			
			return total;
		}
	
	
		// cast List<Double> into double[]
		// http://stackoverflow.com/questions/6018267/how-to-cast-from-listdouble-to-double-in-java
		public static double[] convertDoubleListToArray(List<Double> doubles)
		{		
			double[] target = new double[doubles.size()];
			 for (int i = 0; i < target.length; i++) {
			   // target[i] = doubles.get(i).doubleValue();  // java 1.4 style
			    // or:
			    target[i] = doubles.get(i);                // java 1.5+ style (outboxing)
			 }
			
			
			 return target;
		}
		
		
		public static List<Integer> castDoubleArrayToIntList(double[] doubles)
		{		
			List<Integer> intList = new ArrayList<Integer>();

			for (int i = 0; i < doubles.length; i++) intList.add((int)doubles[i]);
			
			
			return intList;
		}
		
		public static double[] castIntArrayToDouble(int[] intArray)
		{		
			double[] doubleArray = new double[intArray.length];

			for (int i = 0; i < intArray.length; i++) doubleArray[i] = intArray[i];
			
			
			return doubleArray;
		}
		
		public static List<Double> castDoubleArrayToDoubleList(double[] doubles)
		{		
			List<Double> intList = new ArrayList<Double>();

			for (int i = 0; i < doubles.length; i++) intList.add(doubles[i]);
			
			
			return intList;
		}
		
		public static List<byte[]> cast2DByteArrayToByteArrayList(byte[][] arrays)
		{		
			List<byte[]> listArray = new ArrayList<byte[]>();

			for (int i = 0; i < arrays.length; i++) listArray.add(arrays[i]);
			
			
			return listArray;
		}
		
		public static int[] convertIntListToArray(List<Integer> list)
		{

			int[] array = new int[list.size()];
			for(int i = 0; i < list.size(); i++) array[i] = list.get(i);
			
			return array;
		}	
		
		public static byte[] convertByteListToArray(List<Byte> list)
		{

			byte[] array = new byte[list.size()];
			for(int i = 0; i < list.size(); i++) array[i] = list.get(i);
			
			return array;
		}
		
		
		public static byte[] convertHaplotypeListToByteArray(List<Byte> list, boolean constrain01)
		{

			byte[] array = new byte[list.size()];
			byte allele;
			for(int i = 0; i < list.size(); i++) 
			{
				// constrain results to be 0-1, if that was requested
				allele = list.get(i);
				if(constrain01)
				{
					if(allele < 0) allele = 0; // noncalls are assumed to be references
					else if(allele > 1) allele = 1; // all alternates are constrained to be 1
				}
				array[i] = allele;
			}
			
			return array;
		}
		
		public static boolean[] convertHaplotypeListToBooleanArray(List<Byte> list)
		{

			boolean[] array = new boolean[list.size()];
			byte allele;
			boolean alleleBoolean;
			for(int i = 0; i < list.size(); i++) 
			{
				// constrain results to be 0-1
				allele = list.get(i);
				
				
					if(allele <= 0) alleleBoolean = false; // noncalls are assumed to be references
					else /*if(allele >= 1) */alleleBoolean = true; // all alternates are constrained to be 1
				
				array[i] = alleleBoolean;
			}
			
			return array;
		}
// ________________________________________________	
		
		
		
	    
	    /**
	     * Compares a specified column in a 2D array. used for 2D array by column sort for ASCENDING sort
	     *
	     */
	    private static class ColumnComparator_Asc implements Comparator<double[]> {
	    	
	    	private int columnIndex;
	    	public  ColumnComparator_Asc(int columnIndex)
	    	{
	    		this.columnIndex = columnIndex;
	    	}
	        @Override
	        public int compare(double[] row1, double[] row2) {
	            return Double.compare(row1[columnIndex], row2[columnIndex]);
	            // or, before Java 7:
	            // return Double.valueOf(row1[1]).compareTo(Double.valueOf(row2[1]));
	        }    
	    }
	    
	    /**
	     * Compares a specified column in a 2D array. used for 2D array by column sort for DESCENDING sort
	     *
	     */
	    private static class ColumnComparator_Desc implements Comparator<double[]> {
	    	
	    	private int columnIndex;
	    	public  ColumnComparator_Desc(int columnIndex)
	    	{
	    		this.columnIndex = columnIndex;
	    	}
	        @Override
	        public int compare(double[] row1, double[] row2) {
	            return Double.compare( row2[columnIndex],row1[columnIndex]);
	            // or, before Java 7:
	            // return Double.valueOf(row1[1]).compareTo(Double.valueOf(row2[1]));
	        }    
	    }
		
	    
	    
		/**
		 * Sorts a Table(matrix) based on its first row in Descending order
		 * 
		 * @param origTable: the original 2D matrix to be sorted
		 * @return: a new 2D Array (does not modify the old one )
		 */
	    public static double[][] sortTableByRow(double[][] origTable) { return sortTableByRow(origTable,0,true); }
	    
		/**
		 * Sorts a Table(matrix) based on one of its rows
		 * 
		 * @param origTable: the original 2D matrix to be sorted
		 * @param sortRowIndex: the row's index to be sorted by
		 * @param descending: if the ordering should be Descending (true) or Ascending (false)
		 * @return: a new 2D Array (does not modify the old one )
		 */
		public static double[][] sortTableByRow(double[][] origTable, int sortRowIndex, boolean descending) 
		{
			// it is not possible to sort an 2D array by column, if there are duplicate values, as we wont be able to match the indices back to the original array
			
			// SOLUTION: rotate the 2D array to 'by column' -> sort it by column -> rotate it back...
			
			// I. rotate it forward (+90 )
			double[][] rotatedTable = rotateMatrix90Degrees(origTable);
			
			// II. sort matrix by specified column
			if(descending) Arrays.sort(rotatedTable, new ColumnComparator_Desc(origTable.length - sortRowIndex -1) ); // use different comparator to achieve ascending/descending order
			else  Arrays.sort(rotatedTable, new ColumnComparator_Asc(origTable.length - sortRowIndex -1) );
			
			// III. rotate it back -90
			return rotateMatrix_minus_90Degrees(rotatedTable);
		}
		
		
		
		/**
		 * Sorts a Table(matrix) based on its first col in Descending order
		 * 
		 * @param origTable: the original 2D matrix to be sorted
		 * @return: a new 2D Array (does not modify the old one ) 
		 */
	    public static double[][] sortTableByCol(double[][] origTable) { return sortTableByCol(origTable,0,true); }
	    
	    
		/**
		 * Sorts a Table(matrix) based on one of its cols
		 * 
		 * @param origTable: the original 2D matrix to be sorted
		 * @param sortColIndex: the row's index to be sorted by
		 * @param descending: if the ordering should be Descending (true) or Ascending (false)
		 * @return: a new 2D Array (does not modify the old one )
		 */
		public static double[][] sortTableByCol(double[][] origTable, int sortColIndex, boolean descending) 
		{

			// I. because the sort modifies the original array, we need to duplicate it first
	    	double[][] result = new double[origTable.length][];
	    	for(int i = 0; i < origTable.length; i++)	result[i] = origTable[i].clone();
		    	
			// II. sort matrix by specified column
			if(descending) Arrays.sort(result, new ColumnComparator_Desc(sortColIndex) ); // use different comparator to achieve ascending/descending order
			else  Arrays.sort(result, new ColumnComparator_Asc(sortColIndex) );
			

			return result;
		}

	    /**
	     * Rotates a matrix by -90 degrees, so that the LAST COLUMN will become the FIRST ROW
	     * 
	     * @param matrix: the matrix to be rotated
	     * @return a new rotated matrix
	     */
	    public static double[][] rotateMatrix_minus_90Degrees(double[][] matrix)
	    {
	    	// because reverse rows modifies the original array, we need to duplicate it first
	    	double[][] result = new double[matrix.length][];
	    	for(int i = 0; i < matrix.length; i++)	result[i] = matrix[i].clone();
	    	
	    	// to rotate by -90: reverse rows, then transpose
	       	reverseRows( result) ;
	       	result = transposeMatrix(result);
	    	return result;
	    } 

	    
	    /**
	     * Rotates a matrix by +90 degrees, so that the FIRST ROW will become the LAST COLUMN
	     * 
	     * @param matrix: the matrix to be rotated
	     * @return a new rotated matrix
	     */
	    public static double[][] rotateMatrix90Degrees(double[][] matrix)
	    {
	    	// to rotate by +90: transpose, then reverse rows
	    	double[][] result = transposeMatrix(matrix); // creates new array
	    	reverseRows( result) ; // modifies in place
	    	return result;
	    }
	    
	    // http://stackoverflow.com/questions/26197466/transposing-a-matrix-from-a-2d-array
	    public static double[][] transposeMatrix(double[][] matrix)
	    {
	        int m = matrix.length;
	        int n = matrix[0].length;
	        int x;
	        int y;
	        double[][] trasposedMatrix = new double[n][m];

	        for(x = 0; x < n; x++)
	        {
	            for(y = 0; y < m; y++) trasposedMatrix[x][y] = matrix[y][x];
	        }

	        return trasposedMatrix;
	    }
	    
	    
	    public static void reverseRows (double[][] matrix)
	    {
	    	for (int i =0; i < matrix.length; i++) reverseArray(matrix[i]);
	    	
	    }
	    
	    public static int[][] transposeMatrix(int[][] matrix)
	    {
	        int m = matrix.length;
	        int n = matrix[0].length;
	        int x;
	        int y;
	        int[][] trasposedMatrix = new int[n][m];

	        for(x = 0; x < n; x++)
	        {
	            for(y = 0; y < m; y++) trasposedMatrix[x][y] = matrix[y][x];
	        }

	        return trasposedMatrix;
	    }
	    
	    
	    
	 // http://stackoverflow.com/questions/2137755/how-do-i-reverse-an-int-array-in-java
	    public static void reverseArray (double[] someArray)
	    {
	    	double temp;
	    	for(int i = 0; i < someArray.length / 2; i++)
	    	{
	    	    temp = someArray[i];
	    	    someArray[i] = someArray[someArray.length - i - 1];
	    	    someArray[someArray.length - i - 1] = temp;
	    	}
	    }
		
		
		
		
		
		
		
		
		
		/* 
		// this doesnt work... it wasnt possible to circumvent the duplicate element issue...
		public static double[][] sortTableByRow(double[][] origTable, int sortRowIndex, boolean descending) 
		{
			int i;
			int j;
			

			// 0. need to prefill the rows of the new table with NaNs in order to be able to detect gaps
			double[][] newTable = new double[origTable.length][origTable[0].length]; // create new table same size as old one
			for(i = 0; i < newTable.length; i++) Arrays.fill(newTable[i], Double.NaN);

			
			// need to associate each entry with a unique symbol
			// and that unique symbol is then associated with the actual value
			
			
			// I. store the original indices of each entry in the target row
			// Note: for duplicate entries this will only store them once, which will cause 'sparse' arrays later that need to be manually fixed
			LinkedHashMap<Double, Integer> origIndices = new LinkedHashMap<Double, Integer>();
			for(i = 0; i < origTable[sortRowIndex].length; i++) origIndices.put(origTable[sortRowIndex][i], i);
				
		
			
			// II. duplicate target row contents into the new array: http://stackoverflow.com/questions/5785745/make-copy-of-array-java
			System.arraycopy( origTable[sortRowIndex], 0, newTable[sortRowIndex], 0, origTable[sortRowIndex].length );
				//for(i = 0; i < origTable.length; i++)
				//	System.arraycopy( origTable[i], 0, newTable[i], 0, origTable[i].length );	

			
			// III. sort the target row
			Arrays.sort( newTable[sortRowIndex] );
			if(descending)ArrayUtils.reverse(newTable[sortRowIndex] );
			
			
			// IV. Copy over the data from the origTable into the new table at matching indices according to the original table
			int origIndex;
			for(i = 0; i < newTable[sortRowIndex].length; i++) // go through the elements of the new sorted row
			{
				
				origIndex = origIndices.get(newTable[sortRowIndex][i]); // get the original index for this particular entry before sorting
				
				for(j = 0; j < newTable.length; j++) // go through the other rows of the table
				{
					if(sortRowIndex == j) continue; // if its the target row then dont do anything, as that is what we wanted to sort by...
					
					// write in the new table the data from the matching index from the old table
					newTable[j][i] = origTable[j][origIndex];
				}
			
			}
			
			// V. For duplicate entries the above solution will only fill in the FIRST occurence with the original value, the remaining will be left as NaNs
			
			for(i = 0; i < newTable.length; i++) // go through the elements of the new sorted row
			{
				if(i == sortRowIndex) continue; // if its the target row then dont do anything
				
				for(j = 0; j < newTable[i].length; j++) // go through the rows of the orig table
				{	
					if(  Double.isNaN(newTable[i][j]) ) newTable[i][j] = newTable[i][j-1];
				}
			
			}
			
			
			return newTable;
			
		}*/
		
		//  what happens if 2 entries are the same??
		// [10.1, 11.3, 12.5, 12.5, 16.7]
		// [1   ,    2,    3,    3,    4]
		// in the Dic both will be stored as dic[12.5] = 3
		// so when th sorted array will be generated:
		// [10.1, 11.3, 12.5,  NaN, 16.7]
		// there will be left gaps
		
		// Idea ? fill gaps by the 
		
}
