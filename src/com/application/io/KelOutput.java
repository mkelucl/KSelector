package com.application.io;

import java.util.ArrayList;
import java.util.List;

public class KelOutput {
	
	
	// X: [variants][samples]
	public static List<String> produceMapOutput(double[][] X, double[] Y) 
	{
		int i = 0;

		List<String> lines = new ArrayList<String>();
		StringBuilder sb;
		
		// Add header, this should look like this; 	SNP1	SNP2	SNP3	SNP4	....	Y
		sb = new StringBuilder();
		for(i = 0; i < X.length; i++) sb.append("SNP" + (i+1) + "\t" );	
		sb.append("Y");
		lines.add(sb.toString() );
		
		int j;
		
		// add data
		for(i = 0; i < X[0].length; i++) // loop by column...
		{
			sb = new StringBuilder();
			for(j = 0; j < X.length; j++)
			{
				sb.append(X[j][i] + "\t" );
			}
			
			sb.append(Y[i]); // finally add the outcome as the last column
			lines.add(sb.toString() );
		
		}
		
		
		
		return lines;
	}
}
