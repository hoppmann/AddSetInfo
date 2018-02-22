package de.hoppmann.addSetInfo;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;


public class AddSetInfo {
	///////////////////////////
	//////// variables ////////
	///////////////////////////
	
	
	//// logger
	private Log log = Log.getInstance(AddSetInfo.class);
	
	//// variant dependent variables
	private AbstractFeatureReader<VariantContext, LineIterator> reader;
	private VariantContextWriter writer;
	private VCFHeader header;
	
	
	//// naming variables
	private String setField = "set";
	private String consensus = "NCALL";
	private String intersect = "Intersection";
	private String caller = "CALLER";
	
	
	//// input variables
	private int maxN;
	private int minN;
	private String fileInPath;
	private String fileOutPath;
	

	
	
	/////////////////////////////
	//////// constructor ////////
	/////////////////////////////

	
	
	
	/////////////////////////
	//////// methods ////////
	/////////////////////////
	public static void main(String[] args) throws Exception {
		

		AddSetInfo add = new AddSetInfo();
		add.getArguments(args);
		add.readVcf();
		add.prepareWriter();
		add.processVariants();
		add.cleanUp();
		
		
		
	}
	

	
	
	
	
	
	
	
	
	/////////////////
	//////// clean up
	
	private void cleanUp() throws Exception {
		reader.close();
		writer.close();
		log.info("Set values successfully added");

	}
	
	
	
	
	
	
	/////////////////////////
	//////// process variants
	
	//// add set to info field
	private void processVariants() throws IOException {
		
		
		// make log entry
		log.info("Processing variants");
		
		// run over each variant
		Iterator<VariantContext> i = reader.iterator();
		
		
		while (i.hasNext()) {
			VariantContext variant =  i.next();
			processCurVariant(variant); 
		}
	}
	

	
	
	
	
	
	
	//////// process single variant
	private void processCurVariant(VariantContext variant) {
		
		// check if consensus calling produced set field if so split 
		if (! variant.hasAttribute(setField)) {
			log.error("No set information found.");
			System.exit(1);
		}
		String[] set = variant.getAttributeAsString(setField, null).split("-");
		String setEntry = variant.getAttributeAsString(setField,null);
		

		/*
		 * if multiple caller called file process else scip
		 * necessary due to a bug in GATK
		 * check for number of entries and string intersect,
		 * intersect -> all caller called this variant. 
		 */
		
		if (set.length >= minN  && ! set[0].equals(intersect)) {

			// get names of individuals and iterate variant processing over all individuals
			Object[] sampleNames = variant.getSampleNames().toArray();

			
			/*
			 *  iterate over all patIDs and add caller info to format section
			 *  if is called add set info 
			 */
			List<Genotype> gts = new LinkedList<Genotype>();
			for (Object curName : sampleNames) {
				
				// get genotype info
				Genotype geno = variant.getGenotype((String) curName);
				
				// check if called if so add genotype info
				if (geno.isCalled()){
				geno = new GenotypeBuilder(geno).attribute(consensus, set.length).make();
				geno = new GenotypeBuilder(geno).attribute(caller, setEntry).make();
				}
				
				gts.add(geno);
				
			}
			
			writeVariant(variant, gts);
			
			
			
			
			
			// handle the case that all caller found a variable "intersect"
		} else if (set[0].equals(intersect)) {
			
			
			// get names of individuals
			Object[] sampleNames = variant.getSampleNames().toArray();

			/*
			 * iterate over each patID and add caller info to format section
			 */
			List<Genotype> gts = new LinkedList<Genotype>();
			for (Object curName : sampleNames) {

				// get genotype info
				Genotype geno = variant.getGenotype((String) curName);
				
				// check if pat was called if so add total number of callers
				if (geno.isCalled()) {
					geno = new GenotypeBuilder(geno).attribute(consensus, maxN).make();
					geno = new GenotypeBuilder(geno).attribute(caller, setEntry).make();
				}
				
				// save in list for later writing
				gts.add(geno);
			}

			// write genotypes
			writeVariant(variant, gts);

		}
	}

	
	
	
	
	
	
	
	
	///////////////////////////
	//////// write new VCF-File
	
	private void prepareWriter () {

		// add header line for new Set info
		header.addMetaDataLine(new VCFFormatHeaderLine(consensus, 1, VCFHeaderLineType.Character, "Number of callers called this variable"));
		header.addMetaDataLine(new VCFFormatHeaderLine(caller, 1, VCFHeaderLineType.Character, "Name list of caller that called the variant"));
		
		// don't index as this needs dictionary
		VariantContextWriterBuilder builder  = new VariantContextWriterBuilder().unsetOption(Options.INDEX_ON_THE_FLY);
		writer = builder.setOutputFile(fileOutPath).build();
		
		// write out header
		writer.writeHeader(header);
	}

	
	
	private void writeVariant(VariantContext variant, List<Genotype>gts) {
		
		VariantContextBuilder variantBuilder =  new VariantContextBuilder(variant);
		variantBuilder.genotypes(gts);
		VariantContext var2 = variantBuilder.make();
		writer.add(var2);
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/////////////////////////
	//////// read in vcf file
	
	public void readVcf () {
		
		File fileIn = new File(fileInPath);
		if (fileIn.exists()) {
			log.info("File " + fileIn + " found!");
		} else {
			log.error("File " + fileIn + " not found");
			System.exit(1);
		}
			
			//// read in VCF file
			reader = AbstractFeatureReader.getFeatureReader(fileIn.getAbsolutePath(), new VCFCodec(), false); 
			header = (VCFHeader) reader.getHeader();
			
	}
	
	
	
	
	
	
	//////////////////////////
	//////// read in arguments
	
	public void getArguments(String[] args) {
		
		
		if (args.length != 4) {
			log.error("Usage: AddSetInfo <vcf-in> <vcf-out> <N caller used> <min N caller for consensus>");
			System.exit(1);
		}
		
		// get file path
		fileInPath = args[0];
		
		// get file out name
		fileOutPath = args[1];
		
		
		// get max number of caller
		maxN = Integer.parseInt(args[2]);
		
		
		// get the number of caller needed for a call
		minN = Integer.parseInt(args[3]);
		
		
		
	}
	
	
	
	
	
	
	/////////////////////////////////
	//////// getter / setter ////////
	/////////////////////////////////
	
	
	
	
}
