package com.application.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import com.application.utils.HelperMethods;

public class SNPh2Parser {
	/**
	 * @param inputFileLocation: this is a headerless file with following signature: SNP1	0.017952999999999997  ... where each individual SNPs h2 contribution is the 2nd col
	 * @param variantNames: optional parameter, that will filter SNPs by their ID (IE if there are less Variants for some reacion (maybe QC) than what we originally got
	 * @return: a list of the variants LociIndices (IE 0 based lookups)
	 */
	public static double[] parseSNPh2Contr(String inputFileLocation, List<String> variantNames)
	{

		Scanner s;
		String nextItem;
		String[] processedLine;
		
		boolean firstLine = false;
		List<Double> h2Values = new ArrayList<Double>();
		try {
			s = new Scanner (new BufferedReader( new FileReader(inputFileLocation)));
			//s.useDelimiter(System.getProperty("line.separator")); // OS independent line delimiter
			while(s.hasNext() ) // keep going through the file
			{
				nextItem = s.nextLine(); // SNP1	0.017952999999999997	
				processedLine = nextItem.split("\t");
				
				// it is possible that not all SNPs are used, so if such a list was supplied, and this current h2 contribution's SNP is not on the list, then skip it
				if(variantNames != null && variantNames.indexOf(processedLine[0]) == -1)  continue;
				h2Values.add(HelperMethods.tryParseDouble(processedLine[1]));

			}
		} catch (FileNotFoundException e) {System.out.println(inputFileLocation); e.printStackTrace(); }
		
		// convert List to Array in 1 line: http://stackoverflow.com/questions/6018267/how-to-cast-from-listdouble-to-double-in-java
		return h2Values.stream().mapToDouble(d -> d).toArray();  //identity function, Java unboxes automatically to get the double value
	
	}
}
