package com;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;






























import java.util.Scanner;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.FastMath;
import org.ojalgo.access.Access2D.Factory;
import org.ojalgo.matrix.BasicMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;
























import com.application.io.KelOutput;
import com.application.io.SNPh2Parser;
import com.application.io.TextFileWriter;
import com.application.logic.LASSO;
import com.application.logic.OLS;
import com.application.logic.PerformanceSummary;
import com.application.utils.DebugUtil;
import com.io.ldak.LDAKIO;
import com.io.plink.PLINKIO;
import com.io.plink.PLINKIO_BIN;
import com.utils.HelperMethods;


public class KSelectorApp {
	// Master variable to control Debug
	public static Boolean DEBUG = false; // master switch to enable debugging 
	public static Boolean DEBUG_SPEED = false; // if speed should be tested
	//public static Boolean DEBUG_CONSOLE = false; // if the console should be enabled

	// Commands from commandline
	Commands commandlineparser = new Commands();
	private LinkedHashMap<String, Object> parameters;
	
	
	
		// Debug
		private DebugUtil debugger;

		
		public KSelectorApp()
		{
			debugger = new DebugUtil();
		}
		
		
		
		/**
		 * Inits app, parses parameters from command line
		 * and depending on instructions the app is executed
		 */
		public void init(String[] args)
		{ if(DEBUG_SPEED) DebugUtil.instance.initDebugLogic();
		
			// get commands
			parameters = commandlineparser.getParameters(args);
			
			// stop if there are errors from the command line
			if(parameters.containsKey(CommandLineParser.ERROR) ) return;
	
		
			// respect unix writer format
			if(parameters.containsKey(Commands.UNIX_FORCED)) TextFileWriter.FORCE_UNIX_FORMAT = true;
			
		

			List<String> variantNames = null;

			// load custom variant names if present
			if(parameters.containsKey(Commands.VAR_NAMES) )	variantNames = PLINKIO.parseMapFile( commandlineparser.getString(Commands.VAR_NAMES,0)); 
			
			// load custom variant names with A2 alleles if present
			List<String> A2Alleles = null;
			if(parameters.containsKey(Commands.PLINK_FORCEREF_ALLELE) )	
			{ System.out.println("loading in reference alleles from .bim file");
				List<List<String>> bimResults = PLINKIO.parseBimFile( commandlineparser.getString(Commands.PLINK_FORCEREF_ALLELE,0)); 
				variantNames = bimResults.get(0);
				A2Alleles = bimResults.get(1);
			}
			
			
			
			
			// Load in any SNP h2 contributions
			double[] snp_priors = null;
			double[] snph2contr = null;
			if(parameters.containsKey(Commands.SNPH2I))
			{
				snph2contr = SNPh2Parser.parseSNPh2Contr(commandlineparser.getString(Commands.SNPH2I,0),variantNames);
				//System.out.println(Arrays.toString(snph2contr));
				
				System.out.println("loaded heritability contribution for number of SNPs: " + snph2contr.length);
				
				if(parameters.containsKey(Commands.SNPH2I_LASSO))
				{
					snp_priors =  SNPh2Parser.parseSNPh2Contr(commandlineparser.getString(Commands.SNPH2I_LASSO,0),variantNames);
				}
				System.out.println("loaded LASSO Priors for each SNP: " + snph2contr.length);		
			}

				
			// load total h2
			double h2 = -1;
			if(parameters.containsKey(Commands.H2) )h2 =  commandlineparser.getDouble(Commands.H2,0);
			System.out.println("input h2: " + h2);
				
			
			// load in any real causals list
			if(parameters.containsKey(Commands.CAUSALS) ) 
			{ System.out.println("Loading in Causals");
	
				// infer if the file is a KSim causals file or an LDAK file (LDAK causals file is usually called "something.effects")
				boolean isKSim =  commandlineparser.getString(Commands.CAUSALS,0).contains(".effects") ? false : true;
			
				causalIndices = LDAKIO.loadCausalVariantList(commandlineparser.getString(Commands.CAUSALS,0),variantNames, isKSim);
				trueCausalIDs = LDAKIO.loadCausalVariantIDs(commandlineparser.getString(Commands.CAUSALS,0), isKSim);
				
				
				System.out.println("casuals loaded");
				
				// if a results list was supplied
				if(parameters.containsKey(Commands.LASSO_RESULTS) && variantNames != null && snph2contr != null) 
				{
					System.out.println("loading in LASSO results");
					List<String> lassoResults = LASSO.load_LASSO_RESULTS(commandlineparser.getString(Commands.LASSO_RESULTS,0));
					
					List<String> lines = PerformanceSummary.compareAgainstSolutions(lassoResults, trueCausalIDs, variantNames, snph2contr,h2);
					
					String outputFileLocation = commandlineparser.getString(Commands.LASSO_RESULTS,0) + ".perf";
					Path outputFilePath = Paths.get(outputFileLocation);
					TextFileWriter.writeFileOutput(lines, outputFilePath);	
					System.out.println("written results into lassoPerformance.txt");
				}
			}

			// MISC1. parse the determined h2 estimate out from a .reml file that LDAK has produced, which is then written into a specified file that can be read into bash
			if( parameters.containsKey(Commands.PARSE_LDAK_H2) ) 
			{ System.out.println("parsing out h2 from LDAK output");
				
				List<Double> h2_list = parseH2FromLDAK(commandlineparser.getString(Commands.PARSE_LDAK_H2,0));
				
				List<String> lines = new ArrayList<String>();
				String h2Name;
				for(int i = 0; i < h2_list.size(); i++)
				{
					h2Name = "H2_"+(i+1);
					if(i == h2_list.size()-1) h2Name = "H2_ALL";
					lines.add(h2Name + "=\"" + h2_list.get(i) + "\";");
				}
				
				// finally add the 2h that is found only in the non background regions
				double[] regionh2s = h2_list.subList(1, h2_list.size()-1).stream().mapToDouble(Double::doubleValue).toArray();
				double nonBGRegionH2 = 0;
				for(int i = 0; i < regionh2s.length; i++) nonBGRegionH2+= regionh2s[i];
					
				lines.add("H2_REGIONS=\"" + nonBGRegionH2 + "\";");
				
				String outputFileLocation = commandlineparser.getString(Commands.PARSE_LDAK_H2,1);
				Path outputFilePath = Paths.get(outputFileLocation);
				TextFileWriter.writeFileOutput(lines, outputFilePath);	
			}
			
			if(parameters.containsKey(Commands.GET_REGION_SNPS))
			{
				System.out.println("getting all non background SNPs from LDAK's multiblup");
				String locationStem = commandlineparser.getString(Commands.GET_REGION_SNPS,0);
				
				getAllBLUPSNPs(locationStem,true );
				
			}

			if(parameters.containsKey(Commands.TRUE_THINNER))
			{
				System.out.println("extracting true causal snps from a total");
				String truecausalLoc = commandlineparser.getString(Commands.TRUE_THINNER,0);
				String regionLoc = commandlineparser.getString(Commands.TRUE_THINNER,1);
				
				extractTruesFromRegion(truecausalLoc,regionLoc );
				
			}
			
			
			
			if(parameters.containsKey(Commands.INPUT_PHENOTYPE) && parameters.containsKey(Commands.INPUT_PLINK) && parameters.containsKey(Commands.H2) )
			{ 
				
				String inFile = commandlineparser.getString(Commands.INPUT_PLINK,0);
			
				int[] regionBoundaries = null;
				double[] regionh2s = null;
				String regionsLocationStem = null;
				String outputname = "lassoresults.csv";
				
				// load in  region bounderies/h2s if any
				if(parameters.containsKey(Commands.LASSO_AUX_INPUTS))
				{	System.out.println("STARTING H2 GUIDED REGIONAL LASSO");
					String remlFileLocation = commandlineparser.getString(Commands.LASSO_AUX_INPUTS,0);
					List<Double> h2_list = parseH2FromLDAK(remlFileLocation);
					// this will be the same list as above, except the first and the last element, as those are the background region and the total
					regionh2s = h2_list.subList(1, h2_list.size()-1).stream().mapToDouble(Double::doubleValue).toArray();
					
					regionsLocationStem = commandlineparser.getString(Commands.LASSO_AUX_INPUTS,1);
					
					regionBoundaries = getAllBLUPSNPs(regionsLocationStem,false);
					
					// create a backup copy of the original regions
					String directory = "";
					String fileNameStem = "";
					String[] locationArray = regionsLocationStem.split("/");
					if(locationArray.length > 1)  { directory = locationArray[0] + "/"; fileNameStem = locationArray[1]; } 
					else fileNameStem = locationArray[0]; // if there was only 1 element, then there wasn't any directory
					
					File srcDir = new File(directory);
					File destDir = new File(directory.split("/")[0] + "_backup/");
					try {FileUtils.copyDirectory(srcDir, destDir);} catch (IOException e) { System.err.println("CANNOT COPY DIRECTORY"); e.printStackTrace();}
	
				}
				else 
				{
					System.out.println("NO REGIONS SUPPLIED: STARTING H2 GUIDED LASSO");
					outputname = "lassoresults_noregions.csv";
				}
				// else { System.err.println("ABORTING! You must specify region h2s/boundaries!!"); return;}
				
				if(snph2contr == null) { System.err.println("ABORTING! You must specify h2 contributions!!"); return;}
			
				System.out.println("trying to load a BINARY file " +  inFile + ".bed");
	
				//System.out.println("variantNames is null: " + (variantNames == null));
				// I. Load input data:
				// load outcomes (Y)
				double[] Y_phenotypes = PLINKIO.parsePheno(commandlineparser.getString(Commands.INPUT_PHENOTYPE,0))[0];
				
				
				// load predictors (X) [variants][samples]
				//double[][] X_variants = PLINKIO.loadDummyTextData(commandlineparser.getString(Commands.INPUT_PLINK,1));
				
				List<List<String>> sampleNames = null;
				List<List<String>> bimResults = null;
				
				bimResults = PLINKIO.parseBimFile(inFile+ ".bim" ); 
				variantNames = bimResults.get(0);
				A2Alleles = bimResults.get(1);
				int numVariants = bimResults.get(0).size();
				
				sampleNames = PLINKIO.parseSampleNamesFamFile(inFile+ ".fam");
				int numIndividuals = sampleNames.get(0).size();
				
				// load predictors (X) [variants][samples]
				double[][] X_variants = PLINKIO_BIN.readPLINK_Bed(inFile+ ".bed",numIndividuals, numVariants);
			//	LASSO.multiplySNPPrior(X_variants,snp_priors);
				
				System.out.println("finished  loading a BINARY file ");

		
				double shrinkageFactor = -1;
				if(parameters.containsKey(Commands.SHRINKAGE_FACTOR)) shrinkageFactor = commandlineparser.getDouble(Commands.SHRINKAGE_FACTOR,0);
				System.out.println("SHRINKAGE_FACTOR is: "  + shrinkageFactor + " overshrink to squared: " + (shrinkageFactor == -1));
				
				
				// III. Perform LASSO
				double[] blankBetas = new double[X_variants.length];
				
				double[][] LASSOResults = null;
				
				if(parameters.containsKey(Commands.LASSO_AUX_INPUTS)) LASSOResults = LASSO.perform_h2GuidedLASSO_regional(blankBetas,X_variants,Y_phenotypes, snph2contr,regionh2s,regionBoundaries,h2,shrinkageFactor);
				else LASSOResults = LASSO.perform_h2GuidedLASSO_new(blankBetas,X_variants,Y_phenotypes, snph2contr,h2);

				// output results
				LASSO.output_LASSO(LASSOResults, 4,variantNames, outputname,regionsLocationStem,regionBoundaries); // regionsLocationStem "chunks_amb_test/region"
				
			}
	
			
			if(DEBUG_SPEED) System.out.println("finished: " + DebugUtil.instance.SPEEDTEST_RESULTS);
		}

		
	private void extractTruesFromRegion(String truecausalLoc,String regionLoc) 
	{
		
		// I. load truecausals and region snps
		List<String> trueCausals = new ArrayList<String>();
		Scanner s;
		try { s = new Scanner (new BufferedReader( new FileReader(truecausalLoc)));
		while(s.hasNext() ) trueCausals.add(s.nextLine()) ;
		} catch (FileNotFoundException e) {System.out.println(truecausalLoc); e.printStackTrace(); }
	
		System.out.println("loaded number of trueCausals " + trueCausals.size());
		
		List<String> regionVariants = new ArrayList<String>();
		try { s = new Scanner (new BufferedReader( new FileReader(regionLoc)));
		while(s.hasNext() ) regionVariants.add(s.nextLine()) ;
		} catch (FileNotFoundException e) {System.out.println(regionLoc); e.printStackTrace(); }
		
		System.out.println("loaded number of number of variants from region " + regionVariants.size());
		
		
		
		
		// II. go through, true causals and write the ones out that are in the region
		List<String> trueCusalsInRegion = new ArrayList<String>();
		for(int i = 0; i < trueCausals.size(); i++)
		{
			if(  regionVariants.indexOf(trueCausals.get(i)) !=-1 ) trueCusalsInRegion.add(trueCausals.get(i));
		}
		System.out.println("number of true causals located in the region " + trueCusalsInRegion.size());
		
		// III. write these to disc
	
		Path outputFilePath = Paths.get(regionLoc + ".trues");
		TextFileWriter.writeFileOutput(trueCusalsInRegion, outputFilePath);	
		System.out.println("number of true causals located in the region: "+ trueCusalsInRegion.size()+",written to file " +  outputFilePath.toString());
		
		// also check with tagging SNPs!!!
	}


	// plink.lasso
		
		
		public List<String> trueCausalIDs = null;
		public List<Integer> causalIndices = null;
	
		
		
		
		
		
		
		// combines all SNPs from all the regions from an MLBLUP list of regions
		public int[] getAllBLUPSNPs(String locationStem, boolean writeToDisk)
		{
			Scanner s;
			
			// parse the file location, if it has a directory
			String directory = "";
			String fileNameStem = "";
			String[] locationArray = locationStem.split("/");
			if(locationArray.length > 1)  { directory = locationArray[0] + "/"; fileNameStem = locationArray[1]; } 
			else fileNameStem = locationArray[0]; // if there was only 1 element, then there wasn't any directory
			
			// find out how many regions there were in total:
			String filetoLoad = directory+fileNameStem+"_number.txt";
			int numRegions =0;
			try { s = new Scanner (new BufferedReader( new FileReader(filetoLoad)));
				numRegions = HelperMethods.tryParseInt(s.nextLine());
			} catch (FileNotFoundException e) {System.out.println(filetoLoad); e.printStackTrace(); }
			
			System.out.println("loaded from " + filetoLoad + " / numRegions: " + numRegions);
			
			int[] regionBounderies = new int[numRegions];
			
			// go through each region, and load all their SNPs
			List<String> allPredictors = new ArrayList<String>();
			int predictorsInRegion = 0;
			int previousRegionBoundary = 0;
			
			for(int i = 0; i < numRegions ; i++)
			{
				
				predictorsInRegion = 0;
				filetoLoad = directory+fileNameStem+(i+1);
				try { s = new Scanner (new BufferedReader( new FileReader(filetoLoad)));
				
				 while(s.hasNext() ) 
				 {
				 	allPredictors.add(s.nextLine());
				 	predictorsInRegion++;
				 }
				 if(i != 0) previousRegionBoundary = regionBounderies[i -1];
				 regionBounderies[i] = previousRegionBoundary + predictorsInRegion;
				 
				 System.out.println("region ["+(i+1)+"] has boundary at: " +  regionBounderies[i] + " /contains ("+predictorsInRegion+")");
				 
				} catch (FileNotFoundException e) {System.out.println(filetoLoad); e.printStackTrace(); }	
			}

			
			if(writeToDisk)
			{
				String outputFileLocation = "MBLUPSNPs.csv";
				Path outputFilePath = Paths.get(outputFileLocation);
				TextFileWriter.writeFileOutput(allPredictors, outputFilePath);	
				System.out.println("written all("+allPredictors.size()+") multiblup region SNPS into " +  outputFileLocation);
			}
			return regionBounderies;
		}
		
		
		
		
		
		
		public static List<Double> parseH2FromLDAK(String sourceFile)
		{
			String[] nextItemArray;
			String textLine = null;
			Scanner s;
			String nextItem;
			boolean parseStart = false;
			List<Double> h2 = new ArrayList<Double>();
			
			try {
				s = new Scanner (new BufferedReader( new FileReader(sourceFile)));
				//s.useDelimiter(System.getProperty("line.separator")); // OS independent line delimiter
				while(s.hasNext() ) // keep going through the file
				{
					nextItem = s.nextLine(); 
					if(parseStart) // if we have started to parse
					{
						nextItemArray = nextItem.split(" ");
						h2.add(HelperMethods.tryParseDouble(nextItemArray[1]));
					}
					
					// find the line that comes just before where the heritability components for each partition are listed
					if(nextItem.contains("Component Heritability")) parseStart = true; // make note that we can start parsing


				}
			} catch (FileNotFoundException e) {System.out.println(sourceFile); e.printStackTrace(); }
			
			
			return h2;
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		public static void calcZscore(double[][] genotypeMatrix)
		{

			for(int i = 0; i < genotypeMatrix.length; i++) 
			{ 
				if(StatUtils.variance(genotypeMatrix[i]) == 0 ) {  System.err.println("cannot normalise monomorphic variant!, its index: " + i + ". You must exclude these at QC. (use PLINK or similar)"); return;} 

				genotypeMatrix[i] =  StatUtils.normalize(genotypeMatrix[i]);
			}
		}

		
		// test if OLS results are same as PLINK (IE all the matrices etc, are the right way up
		
		public void outputTestData(double[][] X_variants, double[] Y_phenotypes)
		{
			// compare GWAS Top hits against the real solution
			List<String> lines = KelOutput.produceMapOutput(X_variants,Y_phenotypes);
			String outputFileLocation = "data/StudyData.txt";
			Path outputFilePath = Paths.get(outputFileLocation);
			TextFileWriter.writeFileOutput(lines, outputFilePath);	
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		public void ojAlgoTest(double[][] X_variants, double[] Y_phenotypes)
		{
/*			try {
				System.out.println("starting to wait 10 secs");
			    Thread.sleep(10000);                 //1000 milliseconds is one second.
			} catch(InterruptedException ex) {
			    Thread.currentThread().interrupt();
			}
			System.out.println("wait over");*/
			
			
			int numtests = 1000;
			
			final Factory<PrimitiveMatrix> tmpFactory = PrimitiveMatrix.FACTORY; // Initialise ojAlgo
			final BasicMatrix X_Matrix = tmpFactory.columns(X_variants);  // System.out.println("X_Matrix is: " + X_Matrix.toString());
			final BasicMatrix Y_Matrix = tmpFactory.columns(Y_phenotypes);  // System.out.println("Y_Matrix is: " + Y_Matrix.toString());
			
			
			System.out.println("number of Threads used: " + ManagementFactory.getThreadMXBean().getThreadCount() );

			
			DebugUtil.instance.processNext("loading took");
			
			for ( int i = 0; i < numtests; i++)
			{
				final BasicMatrix Xtransposed_Matrix = X_Matrix.transpose();
				final BasicMatrix XtransposedX_Matrix = Xtransposed_Matrix.multiply(X_Matrix);
				final BasicMatrix XtransposedXInverted_Matrix = XtransposedX_Matrix.invert();
				
				final BasicMatrix XtransposedXY_Matrix = Xtransposed_Matrix.multiply(Y_Matrix); //  Xt * Y
				final BasicMatrix ThetaMLE = XtransposedXInverted_Matrix.multiply(XtransposedXY_Matrix); // ( Xt * X )^-1 * Xt * Y
				
				
			}
			DebugUtil.instance.processNext("HYPER THREAD OJ ALGO TOOK"); // ojAlgo uses all cores, so it is around 3.5 faster when multi threaded
		}
		
		
		
		
		
		
		
		public void normalize2DArrayByColumns(double[][] X_variants)
		{ 
			if(X_variants.length == 0) return;
			
			double[] tempArray = new double[X_variants.length]; // create a temp array that will all samples for 1 variant

			double mean = 0;
			double SD = 0;
			
			for (int j =0; j < X_variants[0].length;j++) // go through columns, 0 to last
			{
				for (int i =0; i < X_variants.length;i++) // go through rows
					tempArray[i] = X_variants[i][j]; // no need to recreate tempArray, as we just overwrite it	

				// get mean and SD of newly constructed array
				SD = FastMath.sqrt(StatUtils.variance(tempArray));
				mean = StatUtils.mean(tempArray);
				
				for (int i =0; i < X_variants.length;i++) // go through it again and normalize
					 X_variants[i][j] = ( X_variants[i][j] - mean) / SD;

			}
		}
		
		public void compareGenotypes (double[][] genotypes1, double[][] genotypes2, List<String> variantNames)
		{System.out.println("compareGenotypes");
			
			if(genotypes1.length != genotypes2.length) System.out.println("Genotypes contain different number of variants!!: " + genotypes1.length + " / vs " + genotypes2.length);
			for(int i = 0; i < genotypes1.length; i++)
			{
				for(int j = 0; j < genotypes1[i].length; j++)
				{
					if(genotypes1[i][j] != genotypes2[i][j]) System.out.println("Genotypes different at variant index: " +i+ "("+variantNames.get(i)+"), values are: " +genotypes1[i][j] + " vs " + genotypes2[i][j] );
				}
				
			}
		}
		
		
}
