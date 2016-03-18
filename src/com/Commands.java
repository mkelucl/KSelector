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


package com;

import java.util.ArrayList;

public class Commands extends CommandLineParser {
	public static String INPUT_VCF = "*inputvcf"; // eg: *inputvcf ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz
	public static String REGION = "*region"; // eg: *-region Y 2655180 2658789
	public static String INPUT_PLINK = "*inputplink"; // eg: *inputplink *inputplink short_qc_v2_train
	public static String INPUT_PHENOTYPE = "*inputphenotype"; // eg: *inputphenotype phenotype.pheno  or  *inputphenotype phenotype.phe
	
	public static String ONLY_KEEP_SNPS = "*onlysnps"; // eg: *onlysnps
	public static String NUM_MISSING_GENOTYPES = "*misgenotypes"; // eg: *misgenotypes 1
	public static String MIN_ALT_CALLS = "*minaltcalls"; // eg: *minaltcalls 5
	
	public static String UNIX_FORCED = "*forceunix";  // eg: -forceunix // IE if when writing output files we force a unix linebreak format
	public static String RANDOM_SEED = "*seed";  // eg: *seed 1 
	

	public static String ANALYSIS = "*analysis"; // eg: *analysis
	
	public static String SNPH2I = "*snp2h2"; // eg: *snp2h2 snph2i.csv
	
	public static String SNPH2I_LASSO = "*snp2h2_las"; // eg: *snp2h_las snph2i.csv
	
	public static String H2 = "*h2"; // eg: *h2 0.5
	public static String VAR_NAMES = "*varnames";  // eg: *varnames inputText.map
	
	
	public static String CAUSALS = "*simcausals"; // eg: *simcausals phenotypes.effects  // *causals will be translated as 50causals ??? WTF
	
	public static String LASSO_RESULTS = "*lassoresults"; // eg: *lassoresults plink.lasso  
	public static String PLINK_FORCEREF_ALLELE = "*forcerefallele";  // eg: *forceallele plinkfile.bim
	
	public static String GET_REGION_SNPS = "*getregsnps";  // eg: *getregsnps chunks/region
	public static String PARSE_LDAK_H2 = "*parseh2"; // eg: *parseh2 heritability.reml h2_1.txt

	public static String LASSO_AUX_INPUTS = "*h2helper"; // eg: *h2helper heritability.reml chunks/region
	
	public static String TRUE_THINNER = "*truethin"; // eg: *truethin truecasuals.csv MBLUPSNPs1.csv
	
	public static String SHRINKAGE_FACTOR = "*shfact"; // eg: *shfact 0.5
	
	@Override
	public void setupCommandparams()
	{
		super.setupCommandparams();
		
		noParamsErrorMsg = "usage: *inputplink myData.map myData.ped *inputphenotype phenotype.pheno (*analysis) (*seed 1) (*forceunix)";
		
		
		//            key                               if required    how many parameters are required(0 for any)      error message if required, or if not enough fields specified
		commands.put(REGION, new ArrayList<Object>() {{ add(false);                  add(3);                        add("optional param region must bespecified as as: *region chrom startPs endPos"); }} );
		commands.put(INPUT_VCF, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify a target vcf as: *inputvcf yourVCF.vcf.gz"); }} );	

		
		commands.put(INPUT_PLINK, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify the stem name of the 3 PLINK binary files (.bed,.bim,.fam) as: *inputplink short_qc_v2_train"); }} );	
		commands.put(INPUT_PHENOTYPE, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify a target phenotype file as: *inputphenotype phenotype.pheno  or  *inputphenotype phenotype.phe"); }} );	

		commands.put(ONLY_KEEP_SNPS, new ArrayList<Object>() {{ add(false);                  add(0);                        add(""); }} );	
		commands.put(NUM_MISSING_GENOTYPES, new ArrayList<Object>() {{ add(false);                  add(1);                        add("optional parameter random seed must be specified as: *misgenotypes 1"); }} );	
		commands.put(MIN_ALT_CALLS, new ArrayList<Object>() {{ add(false);                  add(1);                        add("optional parameter random seed must be specified as: *minaltcalls 5"); }} );	
		
		
		commands.put(UNIX_FORCED, new ArrayList<Object>() {{ add(false);                  add(0);                        add(""); }} );
		commands.put(RANDOM_SEED, new ArrayList<Object>() {{ add(false);                  add(1);                        add("optional parameter random seed must be specified as: *seed 1"); }} );
		
		
		commands.put(ANALYSIS, new ArrayList<Object>() {{ add(false);                  add(0);                        add(""); }} );

		commands.put(SNPH2I, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify individual SNPh2 contributions: *snp2h2 snph2i.csv"); }} );
		
		commands.put(SNPH2I_LASSO, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify individual SNPh2 priors for LASSO as: *snp2h2_las SNP_h2s_pos.csv"); }} );
		
		
		commands.put(H2, new ArrayList<Object>() {{ add(false);                  add(1);                        add("must specify true overall h2 as: *h2 0.5"); }} );
		
		commands.put(VAR_NAMES, new ArrayList<Object>() {{ add(false);                  add(1);                        add("to load in custom variant names, a PLINK text map file must be specified as: *varnames inputText.map"); }} );
		commands.put(CAUSALS, new ArrayList<Object>() {{ add(false);                  add(1);                        add("optional parameter load causals file must be specified as: *simcausals phenotypes.effects"); }} );

		commands.put(LASSO_RESULTS, new ArrayList<Object>() {{ add(false);                  add(1);                        add("optional parameter load lasso results file must be specified as: *lassoresults plink.lasso  "); }} );

		commands.put(PLINK_FORCEREF_ALLELE, new ArrayList<Object>() {{ add(false);                  add(1);                        add("to load in custom variant names along with the A2 alleles, a PLINK .bim file must be specified as: *forceallele plinkfile.bim"); }} );
		
		commands.put(GET_REGION_SNPS, new ArrayList<Object>() {{ add(false);                  add(1);                        add("to produce a total list of all SNPs in a BLUP, the regions must be specified as: *getregsnps chunks/region"); }} );
		commands.put(PARSE_LDAK_H2, new ArrayList<Object>() {{ add(false);                  add(2);                        add("to parse out the h2 value from LDAK, input and output must be specified as: *parseh2 heritability.reml h2_1.txt"); }} );	
		
		commands.put(LASSO_AUX_INPUTS, new ArrayList<Object>() {{ add(false);                  add(2);                        add("to specify LDAK region h2s and their variants, these must be specified as: *h2helper heritability.reml chunks/region"); }} );	
		
		commands.put(TRUE_THINNER, new ArrayList<Object>() {{ add(false);                  add(2);                        add("to extract list of true causals existing within a list of snps, these must be specified as: *truethin truecasuals.csv MBLUPSNPs1.csv"); }} );	
		
		
		commands.put(SHRINKAGE_FACTOR, new ArrayList<Object>() {{ add(false);                  add(1);                        add("to set a shrinkage factor, it must be specified as: *shfact 0.5"); }} );
		
		
		
		
		// Handle if Test variant is a LIST or single object
		
		//commands.put(REGION, new ArrayList<Object>() {{  add(false); add(2);  add("String"); add("region must bespecified as as: -region chrom startPs endPos"); }} );
	}
	
	
	
	
}
