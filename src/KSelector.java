//Copyright (c) 2015, Marton Kelemen
//All rights reserved.

//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the 
//documentation and/or other materials provided with the distribution.

//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
//CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
//PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
//WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
//ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import java.awt.EventQueue;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.*;


public class KSelector {

	// Application
	private int _version = 1;   // the version of the client Application. If this doesnt match to what the gameAuth SWF says it wont run
	private Boolean _allowedToRun = false;
	
	private KSelectorApp _application;
	
	
	public static void main(String[] args) {
		
		EventQueue.invokeLater(new Runnable() {
			public void run() { 
				try {
					System.out.println("RECEIVED COMMANDS:");
					for(int i = 0; i < args.length; i++) System.out.println(args[i]);
					System.out.println("__________________");
					
					KSelector window = new KSelector( args);

					
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
		
	}
	


	
	/**
	 * Create the application.
	 */
	public KSelector(String[] args) {
		//System.out.println("main");   
		_application = new KSelectorApp();
		
		if(KSelectorApp.DEBUG)
		{
			// h2 guided LASSO
		//	String[] debugArguments = {"*forcerefallele", "data/ldaksim.bim", "*inputplink", "data/ldaksimTEXT.map", "data/ldaksimTEXT.ped", "*inputphenotype", "data/phenotypes.pheno", "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.48"};  //(from-to "50744624", "50763781")

		//	String[] debugArguments = {"*forcerefallele", "data/keep1000_qc.bim", "*inputplink", "data/keep1000_qcTEXT.map", "data/keep1000_qcTEXT.ped", "*inputphenotype", "data/phenotypes1000.pheno", "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")

		//	String[] debugArguments = {"*forcerefallele", "data/ldaksim.bim", "*inputplink", "data/keep1500TEXT.map", "data/keep1500TEXT.ped", "*inputphenotype", "data/phenotypes1500.pheno", "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")

			
			
			// PLINK lasso analysis
		//	String[] debugArguments = {"*varnames", "data/ldaksimTEXT.map", "*simcausals", "data/phenotypes.effects", "*lassoresults", "data/plink1000.lasso" , "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")

		//	// KSelect lasso analysis
		//	 String[] debugArguments = {"*varnames", "data/keep1000_qcTEXT.map", "*simcausals", "data/phenotypes.effects", "*lassoresults", "data/LASSOresults1000_pruned.txt" , "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")

		//	String[] debugArguments = {"*varnames", "data/keep1500TEXT.map", "*simcausals", "data/phenotypes.effects", "*lassoresults", "data/LASSOresults1500_pruned.txt" , "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")

			// String[] debugArguments = {"*varnames", "data/ldaksimTEXT.map", "*simcausals", "data/phenotypes.effects", "*lassoresults", "data/LASSOresults3000_red.txt" , "*snp2h2", "data/LOO_SNP2_contr.txt" , "*h2" , "0.513406"};  //(from-to "50744624", "50763781")


			/// LASSO KIO
			//String[] debugArguments = {"*inputplink", "data/lasso_small/short_half_test", "*inputphenotype", "data/lasso_small/ksimphen_test.pheno", "*snp2h2", "data/lasso_small/SNP_h2s.csv", "*snp2h2_las", "data/lasso_small/SNP_h2s_pos.csv" , "*h2" , "0.5182"};  //(from-to "50744624", "50763781")

			// LASSO helper stuff
			// String[] debugArguments = {"*h2helper", "amblup.reml", "chunks_amb/region"}; 
			 
			// LASSO: parse h2 in non background regions
			//String[] debugArguments = {"*parseh2", "amblup.reml", "h2.txt"}; 

			
			// select true causals from region
			//String[] debugArguments = {"*truethin", "truecasuals.csv", "lassoresults.csv"};
			
			
			// regional test
			String[] debugArguments = {"*inputplink", "short_train_mblupregions", "*inputphenotype", "phenotypes5_train.pheno", "*snp2h2", "SNPH2s.csv", "*h2helper", "amblup.reml", "chunks_amb/region", "*h2" , "0.210547", "*shfact", "0.5"}; 
			
			
			
			
			_application.init(debugArguments);
		}
		
		else _application.init(args);
		

		
		
	}

}
