package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JToggleButton;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.MathArrays;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.io.RawSequences;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.GammaProposal;
import statalign.mcmc.GaussianProposal;
import statalign.mcmc.HyperbolicPrior;
import statalign.mcmc.InverseGammaPrior;
import statalign.mcmc.McmcCombinationMove;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.UniformPrior;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.structalign.AlignmentMove;
import statalign.model.ext.plugins.structalign.ContinuousPositiveStructAlignMove;
import statalign.model.ext.plugins.structalign.Funcs;
import statalign.model.ext.plugins.structalign.HierarchicalContinuousPositiveStructAlignMove;
import statalign.model.ext.plugins.structalign.LibraryMove;
import statalign.model.ext.plugins.structalign.LinearLink;
import statalign.model.ext.plugins.structalign.MultiNormCholesky;
import statalign.model.ext.plugins.structalign.QuadraticLink;
import statalign.model.ext.plugins.structalign.RotationMove;
import statalign.model.ext.plugins.structalign.RotationProposal;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface;
import statalign.model.ext.plugins.structalign.TranslationMove;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.plugins.structalign.RmsdTrace;
import statalign.postprocess.plugins.structalign.StructTrace;
import statalign.postprocess.plugins.structalign.StructTreeVisualizer;
import statalign.utils.LinkFunction;

public class StructAlign extends ModelExtension implements ActionListener {
	
	
	public double[][][] coordDiffs;
	
	/** The command line identifier of this plugin */
	//private static final String CMD_LINE_PLUGIN_ID = "structal";
	private final String pluginID = "structal";	
	
	@Override
	public String getPluginID() {
		return pluginID;
	}
		
	
	JToggleButton myButton;
	
	public boolean globalSigma = true;
	public boolean useLibrary = false;
	public boolean fixedEpsilon = false;
	public boolean fixedSigma2 = false;
	public boolean fixedTranslation = true;
	
	/** 
	 *  If globalSigma = false then this switches on a spike prior at sigma2Hier. 
	 *  This can also be switched on via a command-line option.
	 * */
	public boolean globalSigmaSpike = false; 
	double[] globalSigmaSpikeParams = {1.35,1.1};
	
	public boolean localEpsilon = false;
	
	double structTemp = 1;	
	
	private boolean USE_IN_ALIGNMENT_PROPOSALS = true;
	
	@Override
	public boolean useInAlignmentProposals() {
		return USE_IN_ALIGNMENT_PROPOSALS;
	}
	
	/** Alpha-C atomic coordinate for each sequence and each residue */
	public double[][][] coords;
	
	/** Crystallographic temperature factors, for weighting epsilon. */
	public double[][] bFactors;
		
	/** Alpha-C atomic coordinates under the current set of rotations/translations */
	public double[][][] rotCoords;
	
	/** Axis of rotation for each sequence */
	public double[][] axes;
	/** Rotation angle for each protein along the rotation axis */
	public double[] angles;
	/** Translation vector for each protein */
	public double[][] xlats;
	
	/** The structure used as the reference for rotations. */
	private int refIndex = 0;	

	/** Parameters of structural drift */
	public double[] sigma2;
	public double sigma2Hier;
	public double nu;
	public double tau;
	public double epsilon;
	// TODO Allow starting values to be specified at command line/GUI
	
	/** Pairwise distances implied by current tree topology */
	public double[][] distanceMatrix;
	/** Covariance matrix implied by current tree topology */
	public double[][] fullCovar;
	/** Current alignment between all leaf sequences */
	public String[] curAlign;
	
	public double[][] oldCovar;
	public double[][] oldDist;
	public String[] oldAlign;
	public double oldLogLi;
	
	/** For caching purposes */
	public HashMap<Integer, MultiNormCholesky> multiNorms;	
	public HashMap<Column, MultiNormCholesky> multiNormsLocal;
	private HashMap<Integer, MultiNormCholesky> oldMultiNorms;
	public HashMap<Column, MultiNormCholesky> oldMultiNormsLocal;
	
	/** Relates the structural and sequence evolutionary timescales */
	private LinkFunction<Double> linkFunction; 
	public String linkType = "linear";  
	// TODO change the above public variables to package visible and put 
	// StructAlign.java in statalign.model.ext.plugins.structalign ?
	
	/* Priors */
	private double sigma2PriorShape = 0.001;
	private double sigma2PriorRate = 0.001;
	public PriorDistribution<Double> sigma2Prior;
	boolean sigma2PriorInitialised = false;
	// sigma2Prior will either be InverseGamma or Hyperbolic, depending
	// on whether globalSigma is switched on. It is defined inside the initRun()
	// method.
	private double epsilonPriorShape = 2;//10; //1; //2;
	private double epsilonPriorRate = 2;//50; //5; //2;
	public PriorDistribution<Double> epsilonPrior;
	boolean epsilonPriorInitialised = false;
	
	/** independence rotation proposal distribution */
	public RotationProposal rotProp;
	
	private double tauPriorShape = 0.001;
	private double tauPriorRate = 0.001;
	public InverseGammaPrior tauPrior = new InverseGammaPrior(tauPriorShape,tauPriorRate);

	private double sigma2HPriorShape = 1;
	private double sigma2HPriorRate = 1;
//	public InverseGammaPrior sigma2HPrior = new InverseGammaPrior(sigma2HPriorShape,sigma2HPriorRate);
	HierarchicalContinuousPositiveStructAlignMove sigma2HMove = null;
	public GammaPrior sigma2HPrior = new GammaPrior(sigma2HPriorShape,sigma2HPriorRate);
//	public HyperbolicPrior sigma2HPrior = new HyperbolicPrior();

	private double nuPriorShape = 1;
	private double nuPriorRate = 6; 
	public GammaPrior nuPrior = new GammaPrior(nuPriorShape,nuPriorRate);
	HierarchicalContinuousPositiveStructAlignMove nuMove = null;

	// priors for rotation and translation are uniform
	// so do not need to be included in M-H ratio
	 	
	
	/** Default proposal weights in this order: 
	 *  align, topology, edge, indel param, subst param, modelext param 
	 *  { 35, 20, 25, 15, 10, 0 };
	 */
	private final int pluginProposalWeight = 50;  
	
	//int sigma2Weight = 5; //15;
	int sigma2Weight = 18; // 
	int tauWeight = 10;
	int sigma2HierWeight = 10; // ORIGINAL
	//int sigma2HierWeight = 0;
	int nuWeight = 0; // ORIGINAL
	int nuWeightIncrement = 10; // ORIGINAL
	//int nuWeight = 0;
	//int epsilonWeight = 2;//10;
	int epsilonWeight = 13; //
	int rotationWeight = 2;
	int translationWeight = 2;
	int libraryWeight = 2;
	int alignmentWeight = 2;
	int alignmentWeightIncrement = 0;
	
	/* Weights for combination moves */
	int alignmentRotationWeight = 8; // ORIGINAL
	int alignmentTranslationWeight = 6; // ORIGINAL
	//int alignmentRotationWeight = 15;
	//int alignmentTranslationWeight = 10;
	int alignmentLibraryWeight = 6;
	int sigmaEpsilonWeight = 4; //
	// This is reallocated to sigma2Weight if epsilon is being fixed 
	
	
	/** Starting value for rotation proposal tuning parameter. */
	public final double angleP = 1000;
	/** Starting value for translation proposal tuning parameter. */
	public final double xlatP = .1;
	
	/** Value to fix sigma at if we're not estimating it. */
	public double fixedSigma2Value = 0.0;
	/** Minimum value for epsilon, to prevent numerical errors. */
	public double MIN_EPSILON = 0.01;
	/** Value to fix epsilon at if we're not estimating it. */
	public double fixedEpsilonValue = 0.0;

	// reference to the postprocessing plugin
	private StructTrace structTrace;
	private RmsdTrace rmsdTrace;
	private StructTreeVisualizer structTree;
	public boolean printRmsd = false;
	public boolean showStructTree = false;

	public StructAlign() {
		// By default this plugin is unselectable unless we read in 
		// coordinate data.
		selectable = false;
	}
	
	@Override
	public List<JComponent> getToolBarItems() {
		myButton = new JToggleButton(new ImageIcon(ClassLoader.getSystemResource("icons/protein.png")));
    	myButton.setToolTipText("Structural alignment mode");
    	myButton.addActionListener(this);
    	myButton.setEnabled(false);
    	myButton.setSelected(false);
    	return Arrays.asList((JComponent)myButton);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		setActive(myButton.isSelected());		
		activateAssociatedPlugins();
	}
	
	private static void printMatrix(double[][] ls) {

		for (int i = 0; i < ls.length; i++) {
			for (int j = 0; j < ls[i].length; j++) {
				System.out.print(ls[i][j] + "\t");
			}
			System.out.println();
		}

	}
	
	//Converts an Integer arraylist to an int array
	public static int[] convertIntegers(List<Integer> integers)
	{
	    int[] ret = new int[integers.size()];
	    Iterator<Integer> iterator = integers.iterator();
	    for (int i = 0; i < ret.length; i++)
	    {
	        ret[i] = iterator.next().intValue();
	    }
	    return ret;
	}
	

	@Override
	public String getUsageInfo() {
		StringBuilder usage = new StringBuilder();
		usage.append("___________________________\n\n");
		usage.append("  StructAlign version 1.0\n\n");
		usage.append("^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n");
		usage.append("java -jar statalign.jar -plugin:structal[OPTION1,OPTION2,...]\n");
		usage.append("OPTIONS: \n");
		usage.append("\tprintRmsd=true\t\t(Prints RMSD and sequence identity for each sample.)\n");
		usage.append("\tprintTree=true\t\t(Prints tree with edge lengths replaced by structural diffusivity.)\n");
		usage.append("\tsigma2=X\t\t(Fixes sigma2 at X)\n");
		usage.append("\tepsilon=X\t\t(Fixes epsilon at X)\n");
		usage.append("\tminEpsilon=X\t\t(Sets minimum value for epsilon to X) [default 0.01]\n");
		usage.append("\tlocalEpsilon\t\t(Uses B-factor information [if available] to scale epsilon per site.)\n");
		usage.append("\tlocalSigma\t\t(Allows each branch to have its own sigma parameter)\n");
		usage.append("\tuseLibrary\t\t(Allows rotation library moves to be used)\n");
		usage.append("\tsigmaSpike\t(Activates spike mixture prior for local sigmas)\n");		
		usage.append("\tsigmaSpike={a_b}\t(Specifies parameters for the Beta prior on spike probability) [default (3.1,1.1) ]\n");		
		usage.append("\tsigma2Prior=PRIOR\t(Sets the prior and hyperparameters for sigma2)\n");
		usage.append("\tepsilonPrior=PRIOR\t(Sets the prior and hyperparameters for epsilon)\n");
		usage.append("\tPRIOR can be one of:\n");
		usage.append("\t\thyp\t\tUses a hyperbolic prior (default)\n");
		usage.append("\t\tg{a_b}\t\tUses a Gamma(a,b) prior\n");
		usage.append("\t\tinvg{a_b}\tUses an InverseGamma(a,b) prior\n");
		usage.append("\t\tunif{a_b}\tUses a Uniform(a,b) prior\n");
		usage.append("\tlink=LINK_FUNCTION\tSets link function between sequence and structure time\n");
		usage.append("\tLINK_FUNCTION can be one of:\n");
		usage.append("\t\tlinear (default)\n");
		usage.append("\t\tquadratic\n");
		usage.append("\nNote that the above syntax is designed to work in bash shells. " +
				"Other shells such as csh may require square brackets to be preceded by a backslash.");
												
		return usage.toString();
	}

	@Override
	public void setActive(boolean active) {
		boolean oldActive = this.active;
		if (oldActive != active) {
			super.setActive(active);
			System.out.println("StructAlign plugin is now "+(active?"enabled":"disabled"));
		}
	}
	
	@Override
	public void setParam(String paramName, String paramValue) {
		if (paramName.equals("epsilon")) {
			fixedEpsilon = true;
			fixedEpsilonValue = Double.parseDouble(paramValue);
			addToFilenameExtension("eps_"+fixedEpsilonValue);
			System.out.println("Fixing epsilon to "+fixedEpsilonValue+".");
		}
		else if (paramName.equals("minEpsilon")) {
			MIN_EPSILON = Double.parseDouble(paramValue);
			addToFilenameExtension("minEps_"+MIN_EPSILON);
			System.out.println("Minimum value for epsilon is now "+MIN_EPSILON+".");
		}
		else if (paramName.equals("sigma2")) {
			fixedSigma2 = true;
			fixedSigma2Value = Double.parseDouble(paramValue);
			addToFilenameExtension("sigma2_"+fixedSigma2Value);
			System.out.println("Fixing sigma2 to "+fixedSigma2Value+".");
		}
		else if (paramName.equals("sigmaSpike")) {
			globalSigmaSpike = true;
			globalSigma = false;
			addToFilenameExtension("spike");
			String[] argString = paramValue.split("\\{",2);
			if (argString[1].endsWith("}")) {				
				String [] args = argString[1].substring(0,argString[1].length()-1).split("_",2);
				if (args.length == 2) {
					System.out.println("Using Gamma("+Double.parseDouble(args[0])+","+Double.parseDouble(args[1])+
							") prior for "+paramName+".");
					globalSigmaSpikeParams = new double[] {Double.parseDouble(args[0]),Double.parseDouble(args[1])};
				}
				else {
					throw new IllegalArgumentException(
							"Spike parameters must be specifed in the form\n-plugin:structal[sigmaSpike={a_b}]\n");
				}
			}
			addToFilenameExtension(globalSigmaSpikeParams[0]+"_"+globalSigmaSpikeParams[1]);
			System.out.println("Activating spike prior for local sigmas.");
		}
		else if (paramName.equals("sigma2Prior")) {
			sigma2Prior = setPrior(paramName,paramValue);
			sigma2PriorInitialised = (sigma2Prior != null);
		}
		else if (paramName.equals("epsilonPrior")) {
			epsilonPrior = setPrior(paramName,paramValue);
			epsilonPriorInitialised = (epsilonPrior != null);
		}
		else if (paramName.equals("link")) {
			linkType = paramValue;
			addToFilenameExtension("link_"+paramValue);
			System.out.println("Setting link function to "+paramValue+".");
		}
		else {
			super.setParam(paramName,paramValue);
		}
	}
	@Override
	public void setParam(String paramName, Number paramValue) {
		if (paramName.equals("epsilon")) {
			fixedEpsilon = true;
			fixedEpsilonValue = (Double) paramValue;
			addToFilenameExtension("eps_"+fixedEpsilonValue);
			System.out.println("Fixing epsilon to "+fixedEpsilonValue+".");
		}
		else {
			super.setParam(paramName,paramValue);
		}
	}
	@Override
	public void setParam(String paramName, boolean paramValue) {
		if (paramName.equals("localSigma")) {
			globalSigma = false;
		}
		else if (paramName.equals("localEpsilon")) {
			localEpsilon = true;
			System.out.println("Using B-factor information to scale epsilon.");
		}
		else if (paramName.equals("printRmsd")) {
			printRmsd = true;			
		}
		else if (paramName.equals("printTree")) {
			showStructTree = true;			
		}
		else if (paramName.equals("useLibrary")) {
			useLibrary = true;
		}
		else if (paramName.equals("sigmaSpike")) {
			globalSigmaSpike = true;
			globalSigma = false;
			addToFilenameExtension("spike");
			addToFilenameExtension(globalSigmaSpikeParams[0]+"_"+globalSigmaSpikeParams[1]);
			System.out.println("Activating spike prior for local sigmas.");
		}
		else {
			super.setParam(paramName,paramValue);
		}
	}
	
	private PriorDistribution<Double> setPrior(String paramName, 
						String paramValue) {
		if (paramValue.startsWith("hyp")) {
			addToFilenameExtension(paramName+"_hyp");
			System.out.println("Using hyperbolic prior for "+paramName+".");
			return new HyperbolicPrior();
		}
		else if (paramValue.startsWith("unif{")) {
			String[] argString = paramValue.split("\\{",2);
			if (argString[1].endsWith("}")) {				
				String [] args = argString[1].substring(0,argString[1].length()-1).split("_",2);
				if (args.length == 2) {
					addToFilenameExtension(paramName+"_u_"+args[0]+"_"+args[1]);
					System.out.println("Using Unif("+Double.parseDouble(args[0])+","+Double.parseDouble(args[1])+
							") prior for "+paramName+".");
					return new UniformPrior(Double.parseDouble(args[0]),Double.parseDouble(args[1]));
				}
				else {
					throw new IllegalArgumentException(
							"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=unif{a_b}]\n");
				}
			}
			else {
				throw new IllegalArgumentException(
					"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=unif{a_b}]\n");
			}
		}
		else if (paramValue.startsWith("g{")) {
			String[] argString = paramValue.split("\\{",2);
			if (argString[1].endsWith("}")) {				
				String [] args = argString[1].substring(0,argString[1].length()-1).split("_",2);
				if (args.length == 2) {
					addToFilenameExtension(paramName+"_g_"+args[0]+"_"+args[1]);
					System.out.println("Using Gamma("+Double.parseDouble(args[0])+","+Double.parseDouble(args[1])+
							") prior for "+paramName+".");
					return new GammaPrior(Double.parseDouble(args[0]),Double.parseDouble(args[1]));
				}
				else {
					throw new IllegalArgumentException(
							"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=g{a_b}]\n");
				}
			}
			else {
				throw new IllegalArgumentException(
					"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=g{a_b}]\n");
			}
		}
		else if (paramValue.startsWith("invg{")) {
			String[] argString = paramValue.split("\\{",2);
			if (argString[1].endsWith("}")) {
				String [] args = argString[1].substring(0,argString[1].length()-1).split("_",2);
				if (args.length == 2) {
					addToFilenameExtension(paramName+"_invg_"+args[0]+"_"+args[1]);
					System.out.println("Using InvGamma("+Double.parseDouble(args[0])+","+Double.parseDouble(args[1])+
							") prior for "+paramName+".");
					return new InverseGammaPrior(Double.parseDouble(args[0]),Double.parseDouble(args[1]));
				}
				else {
					throw new IllegalArgumentException(
							"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=invg{a_b}]\n");
				}
			}
			else {
				throw new IllegalArgumentException(
					"Prior parameters must be specifed in the form\n-plugin:structal[sigma2Prior=invg{a_b}]\n");
			}
		}
		else {
			throw new IllegalArgumentException("Unrecognised prior specification "+paramName+"="+paramValue+".");
		}
	}
	@Override
	public void init() {
	}
	
	@Override
	protected void resetData() {
		super.resetData();
		distanceMatrix = null;		
		fullCovar = null;		
		curAlign = null;
		
		oldCovar = null;
		oldDist = null;
		oldAlign = null;
		oldLogLi = Double.NEGATIVE_INFINITY;
	}
//	@Override
//	public ModelExtension reset() {
//		return new StructAlign();
//	}
	@Override
	public void initRun(InputData inputData) throws IllegalArgumentException {
		super.initRun(inputData);
		HashMap<String, Integer> seqMap = new HashMap<String, Integer>();
		int i = 0;
		for(String name : inputData.seqs.getSeqnames())
			seqMap.put(name.toUpperCase(), i++);
		coords = new double[inputData.seqs.size()][][];
		if (localEpsilon) bFactors = new double[inputData.seqs.size()][];
		for(DataType data : inputData.auxData) {
			if(!(data instanceof ProteinSkeletons))
				continue;
			ProteinSkeletons ps = (ProteinSkeletons) data;
			if (localEpsilon && ps.bFactors.size() == 0) {
				throw new RuntimeException("No B-factor data available: cannot use localEpsilon mode.");
			}
			for(i = 0; i < ps.names.size(); i++) {				
				String name = ps.names.get(i).toUpperCase();
				if(!seqMap.containsKey(name))
					throw new IllegalArgumentException("structalign: missing sequence or duplicate structure for "+name);
				int ind = seqMap.get(name);
				int len = inputData.seqs.getSequence(ind).replaceAll("-", "").length();
				List<double[]> cl = ps.coords.get(i);
				List<Double> bF = ps.bFactors.get(i);
				if (localEpsilon && bF.size() == 0) {
					localEpsilon = false;
					System.out.println("No B-factor data available for "+name+": cannot use localEpsilon mode.");					
				}
				if(len != cl.size())
					throw new IllegalArgumentException("structalign: sequence length mismatch with structure file for seq "+name);
				coords[ind] = new double[len][];
				if (localEpsilon) bFactors[ind] = new double[len];
				double bFactorMean = 0;
				for(int j = 0; j < len; j++) {
					 coords[ind][j] = Utils.copyOf(cl.get(j));
					 if (localEpsilon) {
						 bFactors[ind][j] = bF.get(j);
						 bFactorMean += bFactors[ind][j]/len;
					 }
				}
				if (localEpsilon) {
					for(int j = 0; j < len; j++) {
						if (bFactorMean == 0) {
							bFactors[ind][j] = 1;	
						}
						else {
							bFactors[ind][j] /= bFactorMean;	
						}						
					}
				}
				// center all coordinates to mean zero so that rotations are around center of gravity
				RealMatrix temp = new Array2DRowRealMatrix(coords[ind]);
				RealVector mean = Funcs.meanVector(temp);
				for(int j = 0; j < len; j++)
					 coords[ind][j]= temp.getRowVector(j).subtract(mean).toArray();
				seqMap.remove(name);
			}
		}
		while (coords[refIndex]==null) ++refIndex;	
		int nStructures=0; 
		for (int ii=0; ii<coords.length; ii++) {
			if (coords[ii]!=null) ++nStructures;
		}
		if (nStructures < 2) {
			throw new RuntimeException("Cannot run StructAlign with fewer than two structures.");
		}
		
//		if(seqMap.size() > 0)
//			throw new IllegalArgumentException("structalign: missing structure for sequence "+seqMap.keySet().iterator().next());
		
		if (useLibrary) {
			rotProp = new RotationProposal(this);
		}
		
		if (globalSigma && showStructTree) {
			System.out.println("structTree option can only be used in conjunction with localSigma. ");
		}
		
		rotCoords = new double[coords.length][][];
		axes = new double[coords.length][];
		angles = new double[coords.length];
		xlats = new double[coords.length][];

		axes[refIndex] = new double[] { 1, 0, 0 };
		angles[refIndex] = 0;
		xlats[refIndex] = new double[] { 0, 0, 0 };
						
		multiNorms = new HashMap<Integer, MultiNormCholesky>();	
		multiNormsLocal = new HashMap<Column, MultiNormCholesky>();
		oldMultiNorms = new HashMap<Integer, MultiNormCholesky>();
		oldMultiNormsLocal = new HashMap<Column, MultiNormCholesky>();
		
		// MCMC parameters
		sigma2Hier = 1;
		nu = 1;
		tau = 50;
		if (fixedEpsilon) {
			epsilon = fixedEpsilonValue;
			sigma2Weight += sigmaEpsilonWeight;
		}
		else {
			//epsilon = 100;
			epsilon = 50;
			//epsilon = 20;
		}
		
		if (fixedSigma2) {
			globalSigma = true;
		}
		// number of branches in the tree is 2*leaves - 1
		if (globalSigma) {
			sigma2 = new double[1];
		}
		else {
			sigma2 = new double[2*coords.length - 1];
		}
		
		for(i = 0; i < sigma2.length; i++) {
			sigma2[i] = 1;
		}
		if (fixedSigma2) {
			sigma2[0] = fixedSigma2Value;
			epsilonWeight += sigmaEpsilonWeight;
		}
		
		if (!sigma2PriorInitialised && !fixedSigma2) {
			if (globalSigma) {
				if (fixedEpsilon) {
					sigma2Prior = new GammaPrior(2,2);
				}
				else {
					sigma2Prior = new GammaPrior(1,1);
					//sigma2Prior = new HyperbolicPrior();
				}
			}
			else {
				//sigma2Prior = new InverseGammaPrior(sigma2PriorShape,sigma2PriorRate);
				//sigma2Prior = new LinearPrior();
				//sigma2Prior = new UniformPrior();
				sigma2Prior = new GammaPrior(1,1);
			}
			sigma2PriorInitialised = true;
		}
		
		if (!epsilonPriorInitialised && !fixedEpsilon) {
			epsilonPrior = new GammaPrior(epsilonPriorShape,epsilonPriorRate);
			epsilonPriorInitialised = true;
		}
		
		if(linkType.equals("linear")) {
			linkFunction = new LinearLink();
		}
		else if (linkType.equals("quadratic")) {
			linkFunction = new QuadraticLink(); 
		}
		else {
			throw new IllegalArgumentException("Invalid link function selected.");
		}
		
		// for now, don't allow hierarchical sigma unless link function is linear
		if(!linkType.equals("linear")){
			globalSigma = true;
		}
	}
	
	@Override
	protected void initMcmc(InputData inputData) {
		/* Add alignment and rotation/translation moves */		
		RotationMove rotationMove = new RotationMove(this,"rotation"); 
		addMcmcMove(rotationMove,rotationWeight); 
		
		LibraryMove libraryMove = null;
		if (useLibrary) {
			libraryMove = new LibraryMove(this,"library");
			addMcmcMove(libraryMove,libraryWeight);
		}
		
		if (!inputData.pars.fixAlign) {
			AlignmentMove alignmentMove = new AlignmentMove(this,"alignment");
			addMcmcMove(alignmentMove,alignmentWeight,alignmentWeightIncrement); 
			
			/* Combination moves */
			ArrayList<McmcMove> alignmentRotation = new ArrayList<McmcMove>();
			alignmentRotation.add(alignmentMove);
			alignmentRotation.add(rotationMove);
			McmcCombinationMove alignmentRotationMove = 
				new McmcCombinationMove(alignmentRotation);
			addMcmcMove(alignmentRotationMove,alignmentRotationWeight);
			
			if(!fixedTranslation){
				TranslationMove translationMove = new TranslationMove(this,"translation");
				addMcmcMove(translationMove,translationWeight);
				
				ArrayList<McmcMove> alignmentTranslation = new ArrayList<McmcMove>(); 
				alignmentTranslation.add(alignmentMove);
				alignmentTranslation.add(translationMove);
				McmcCombinationMove alignmentTranslationMove = 
					new McmcCombinationMove(alignmentTranslation);
				addMcmcMove(alignmentTranslationMove,alignmentTranslationWeight);
			}
				
			if (useLibrary) { 
				ArrayList<McmcMove> alignmentLibrary = new ArrayList<McmcMove>();
				alignmentLibrary.add(alignmentMove);
				alignmentLibrary.add(libraryMove);
				McmcCombinationMove alignmentLibraryMove = 
					new McmcCombinationMove(alignmentLibrary);
				addMcmcMove(alignmentLibraryMove,alignmentLibraryWeight);
			}
		}
		
		/** Add moves for scalar parameters */
		StructAlignParameterInterface paramInterfaceGenerator = new StructAlignParameterInterface(this); 
		
		ParameterInterface tauInterface = paramInterfaceGenerator.new TauInterface();
		ContinuousPositiveStructAlignMove tauMove = 
			new ContinuousPositiveStructAlignMove(this,tauInterface,tauPrior,new GammaProposal(0.001,0.001),"tau");
		tauMove.moveParams.setPlottable();
		tauMove.moveParams.setPlotSide(1);
		addMcmcMove(tauMove,tauWeight);
		
		ContinuousPositiveStructAlignMove epsilonMove = null;
		if (!fixedEpsilon) {
			ParameterInterface epsilonInterface = paramInterfaceGenerator.new EpsilonInterface();
			epsilonMove = 
				new ContinuousPositiveStructAlignMove(this,epsilonInterface,epsilonPrior,new GaussianProposal(),"eps");
				//new ContinuousPositiveStructAlignMove(this,epsilonInterface,epsilonPrior,new GammaProposal(0.001,0.001),"eps");
			epsilonMove.setMinValue(MIN_EPSILON);
			epsilonMove.moveParams.setPlottable();
			epsilonMove.moveParams.setPlotSide(1);
			addMcmcMove(epsilonMove,epsilonWeight);
		}
		
		if (!fixedSigma2) {			
			if (!globalSigma) {
				ParameterInterface sigma2HInterface = paramInterfaceGenerator.new Sigma2HInterface();
				sigma2HMove = new HierarchicalContinuousPositiveStructAlignMove(this,sigma2HInterface,sigma2HPrior,new GammaProposal(0.001,0.001),"s2_g");
				sigma2HMove.moveParams.setPlottable();
				sigma2HMove.moveParams.setPlotSide(0);
				addMcmcMove(sigma2HMove,sigma2HierWeight); 
				
				ParameterInterface nuInterface = paramInterfaceGenerator.new NuInterface();
				nuMove = new HierarchicalContinuousPositiveStructAlignMove(this,nuInterface,nuPrior,new GammaProposal(0.001,0.001),"nu");
				nuMove.moveParams.setPlottable();
				nuMove.moveParams.setPlotSide(1);
				nuMove.onlySampleIfAtLeastTwoChildrenNotFixed();
				nuMove.setMinValue(0.1);
				if (!globalSigmaSpike) {
					nuWeight += nuWeightIncrement;
					nuWeightIncrement = 0;
				}
				addMcmcMove(nuMove,nuWeight,nuWeightIncrement);
			}
			
			for (int j=0; j<sigma2.length; j++) {
				if (j<=sigma2.length/2 && coords[j]==null) continue; 
				String sigmaName;
				if (sigma2.length == 1) {
					sigmaName = "s2";
				}
				else {
					sigmaName = "s2_"+j;
				}
				ParameterInterface sigma2Interface = paramInterfaceGenerator.new Sigma2Interface(j);
//				ProposalDistribution prop = null;
//				if (globalSigma) prop = new GaussianProposal();
//				else prop = new MultiplicativeProposal(); 
				ContinuousPositiveStructAlignMove m = new ContinuousPositiveStructAlignMove(
															this,sigma2Interface,
															sigma2Prior,new GaussianProposal(),sigmaName);
															//sigma2Prior,new GammaProposal(0.001,0.001),sigmaName);
				if (!globalSigma && j == sigma2.length - 1) {
					continue;
					// i.e. don't add the last one if we have
					// more than one
				}
				if (j<=sigma2.length/2) {
					// plot only for tips
					m.moveParams.setPlottable(); 
					m.moveParams.setPlotSide(0);
				}				
				addMcmcMove(m,sigma2Weight);
				if (!globalSigma) {
					sigma2HMove.addChildMove(m);
					if (globalSigmaSpike) {
						sigma2HMove.setSpikeParams(globalSigmaSpikeParams);
						sigma2HMove.disallowSpikeSelection();
					}					
					m.addParent(sigma2HMove);
					nuMove.addChildMove(m);
					// Don't add nuMove as a parent, because otherwise
					// we'll double count the prior.
				}
				
				if (sigma2.length == 1 && !fixedEpsilon) {
					ArrayList<McmcMove> sigmaEpsilon = new ArrayList<McmcMove>();
					sigmaEpsilon.add(m);
					sigmaEpsilon.add(epsilonMove);
					McmcCombinationMove sigmaEpsilonMove = 
						new McmcCombinationMove(sigmaEpsilon);
					addMcmcMove(sigmaEpsilonMove,sigmaEpsilonWeight);
				}
			}
		}		
	}
	
	@Override
	public void beforeSampling(Tree tree) {
		/* check for protein with no residues aligned to reference protein, resample alignment if any found */
		boolean stop = false;
		while(!stop){
			String[] align = tree.getState().getLeafAlign();	
			String ref = align[refIndex];
			proteinLoop:
			for(int i = 0; i < align.length; i++){
				if (i==refIndex) continue; 				
				String other = align[i];
				int countAligned = 0;
				for(int j = 0; j < align[refIndex].length(); j++)
					countAligned += (ref.charAt(j) != '-' & other.charAt(j) != '-') ? 1 : 0;
				if(countAligned == 0){
					tree.root.selectAndResampleAlignment();
					break proteinLoop;
				} else if (i == align.length - 1) {stop = true;}
			}
		}
		Funcs.initLSRotations(tree,coords,xlats,axes,angles);
		//calcAllRotations();
		computeLogLikeFactor(tree);
	}
	
	@Override
	public void afterFirstHalfBurnin() {	
		if (!globalSigma && globalSigmaSpike) {
			sigma2HMove.allowSpikeSelection();
			zeroAllMoveCounts();
		}
	}
	public void afterBurnin() {
		if(nuMove!=null){
			nuMove.alwaysSample();
			nuMove.setMinValue(0.001);
		}
	}
	
	public double computeLogLikeFactor(Tree tree) {
		String[] align = tree.getState().getLeafAlign();
		checkConsAlign(align); 		
		curAlign = align;				 
		
		double[][] covar = calcFullCovar(tree);
		checkConsCovar(covar); 
		fullCovar = covar;
		
		if(!checkConsRots() && rotCoords[refIndex] == null)
			calcAllRotations();
		
		double logli = calcAllColumnContrib();
		checkConsLogLike(logli); 
		curLogLike = logli;
				
		return curLogLike;
	}
	
	@Override
	public double logLikeFactor(Tree tree) {		

		// Compute log likelihood if not yet computed
		if (curLogLike==0) return computeLogLikeFactor(tree);

		if (Utils.DEBUG) {
			double oldLogLike = curLogLike;
			if (Math.abs(oldLogLike - computeLogLikeFactor(tree)) > 1e-8) {
				throw new RuntimeException("Inconsistency in logLikeFactor: "+
						oldLogLike +" != "+curLogLike);
			}
		}
			
		// If it's non-zero, then we return its current value
		//return computeLogLikeFactor(tree);	
		return curLogLike;
	}
	
	/**
	 * Calculates the structural likelihood contribution of a single alignment column
	 * @param col the column, id of the residue for each sequence (or -1 if gapped in column)
	 * @return the likelihood contribution
	 */		
	public double columnContrib(int[] _col) {
		// count the number of ungapped positions in the column
		int numMatch = 0;
		int[] col = _col.clone();
		for(int i = 0; i < col.length; i++){
			if (col[i]==0) col[i] = -1;
		}
		for(int i = 0; i < col.length; i++){
			if (coords[i]==null) col[i] = -1;	
			if(col[i]!=-1)
				numMatch++;
		}
		if(numMatch == 0) 
			return 0;
		// collect indices of ungapped positions
		int[] notgap = new int[numMatch];
		int columnCode = 0;
		int j = 0;		
		for(int i = 0; i < col.length; i++)  {
			if(col[i]!=-1) {
				notgap[j++] = i;
				columnCode |= (1 << i); 
				/*if this column has exactly the same indel pattern/column code, the alignment has not changed, and 
				 * can reuse the cholesky decompositions - doesnt matter what column you are in, only matters if you have 
				 * same number of amino acids in a column (same number of means and same size of covariance matrix)
				 */
			}						
		}
		
		/*		 
		 * Under localEpsilon mode, the covariance depends on the column,
		 * not just the indel pattern of the column, but we can still
		 * cache the Cholesky decompositions to be re-used for columns
		 * that do not change (since most of the alignment columns do 
		 * not change during an alignment move, this could still yield
		 * a significant speedup).
		 */			
		
		MultiNormCholesky multiNorm = null;
		if (localEpsilon) multiNorm = multiNormsLocal.get(new Column(col));
		else multiNorm = multiNorms.get(columnCode);				
		MultiNormCholesky multiNorm2 = null;
		if (Utils.DEBUG){
			double[][] subCovar = Funcs.getSubMatrix(fullCovar, notgap, notgap);
			// create normal distribution with mean 0 and covariance subCovar
			if (localEpsilon) addLocalEpsilonToDiagonal(subCovar,notgap,col);
			multiNorm2 = new MultiNormCholesky(new double[numMatch], subCovar);
		}
		
		if (multiNorm == null) {
			// extract covariance corresponding to ungapped positions
			double[][] subCovar = Funcs.getSubMatrix(fullCovar, notgap, notgap);
			if (localEpsilon) addLocalEpsilonToDiagonal(subCovar,notgap,col);
			// create normal distribution with mean 0 and covariance subCovar
			multiNorm = new MultiNormCholesky(new double[numMatch], subCovar);
			if (localEpsilon) multiNormsLocal.put(new Column(col), multiNorm);
			else multiNorms.put(columnCode, multiNorm);
		}		
			
		double logli = 0;
		double[] vals = new double[numMatch];
		double aa1;
		double aa2;
				
		// loop over all 3 coordinates
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++) {
				aa1 = rotCoords[notgap[i]][col[notgap[i]]][j];
				aa2 = rotCoords[notgap[i]][col[notgap[i]]-1][j];
				vals[i] =  coordDiffs[notgap[i]][col[notgap[i]]-1][j] = aa1-aa2;
			}
			
			if (Utils.DEBUG  && multiNorm.logDensity(vals) != multiNorm2.logDensity(vals)) {
				System.out.print("col = [");
				for (int k=0; k<col.length; k++) System.out.print(col[k]+",");
				System.out.println("] ("+columnCode+")");
				for (int key : multiNorms.keySet()) {
					if (multiNorms.get(key).getMeans().length==vals.length) {
						System.out.println(key+" "+multiNorms.get(key).logDensity(vals));
					}
				}
				throw new RuntimeException(
						"Inconsistency: "+multiNorm.logDensity(vals)+" != "
						+multiNorm2.logDensity(vals));
			}
			
			logli += multiNorm.logDensity(vals);
		}
		return logli;
	}
	
	/**
	 * Calculates the structural likelihood contribution of a single alignment column
	 * @param col the column, id of the residue for each sequence (or -1 if gapped in column)
	 * @return the likelihood contribution
	 */		
	public double columnContribBondAlign(int[] col, String[] align) {
		
//		System.out.println(Arrays.toString(col));
//		for (String st : align)
//			System.out.println(st);
		
		int bondLen = 0;
		double logli = 0;
//		int[] inds = new int[col.length];		// current char indices
		ArrayList<Integer> bondLens = new ArrayList<Integer>();
			
		for (int seq = 0; seq < col.length; seq++) {
			if(col[seq] == -1){
				bondLens.add(0);
				continue;
			}
			int alignIdx = -1;
			int aaIdx = -1;
			while(aaIdx != col[seq]){
				alignIdx++;
				if (align[seq].charAt(alignIdx) != '-')
					aaIdx++;
			}
//			System.out.println(aaIdx + " aa");
//			System.out.println(alignIdx + " align");
			if (alignIdx == 0){
				bondLens.add(0);
				continue;
			}
			else{
			//If there is an amino acid in the column, increase the bond length count until you find 
			//the preceding amino acid (or the beginning of the sequence)
				bondLen = 1;
				int i = alignIdx - 1;
				while (i != -1 && align[seq].charAt(i) == '-') {
					bondLen++;
					i--;
				}
				//If there is a bond length estimated across the 5' terminal gap, erase it
				//i == 0 can not be used for the "beginning of the sequence" because if there are amino 
				//acids in the first two columns The second amino acid will have a bond length we want to save and i == 0
				if(i == -1){
					bondLens.add(0);
					continue;
				}
				bondLens.add(bondLen);
			}
		}
		//create a set of all the unique bond lens
		Set<Integer> uniqBonds = new HashSet<Integer>();
		uniqBonds.addAll(bondLens);
		//remove 0 bond lengths (correspond to gaps in the current column)
		uniqBonds.remove(0);
		
		int[] bondLensArr = convertIntegers(bondLens);
		
		if(uniqBonds.size() == 0)
			return 0;
		
		ArrayList<Integer> indices = new ArrayList<Integer>();
		
		for(int bond : uniqBonds){
//			System.out.println("bond length:\t"+bond);
			for(int i = 0; i< bondLensArr.length; i++){
				if (bondLensArr[i] == bond)
					indices.add(i);
			}
//			create array with sequences with this bond length
			int[] seqs = convertIntegers(indices);
//			System.out.println("Seqs with bond length:\t"+Arrays.toString(seqs));
			
			//for every seq with a bond length in this column, do a bitshift operation
			int j = 0;
			int columnCode = 0;
//			System.out.println(Arrays.toString(seqs));
			for(int i = 0; i < bondLensArr.length; i++)  {
				if(j < seqs.length && seqs[j]==i){
					columnCode |= (1 << i);
					j++;
				}
					/*integers are 32 bit, meaning there are  2^32 unique codes, 
					is it possible to create duplicate codes for different indel patterns?*/
			}
//			System.out.println("column code:\t"+columnCode);
			indices.clear();
				/*
				 * if this column has exactly the same indel pattern/column
				 * code, can reuse the
				 * cholesky decompositions - doesnt matter what column you are
				 * in, only matters if you have same sequences with amino acids in a
				 * column (same number of means and same subset of covariance
				 * matrix)
				 */

		/*
		 * Under localEpsilon mode, the covariance depends on the column, not
		 * just the indel pattern of the column, but we can still cache the
		 * Cholesky decompositions to be re-used for columns that do not change
		 * (since most of the alignment columns do not change during an
		 * alignment move, this could still yield a significant speedup).
		 */
	
			MultiNormCholesky multiNorm = null;
			if(multiNorms.get(columnCode) != null)
				multiNorm = multiNorms.get(columnCode);
	
			if (multiNorm == null) {
				// extract covariance corresponding to ungapped positions
				double[][] subCovar = Funcs.getSubMatrix(fullCovar, seqs, seqs);
//				System.out.println(subCovar.length +","+ subCovar[0].length);
				// create normal distribution with mean 0 and covariance subCovar
				multiNorm = new MultiNormCholesky(new double[seqs.length], subCovar);
				multiNorms.put(columnCode, multiNorm);
			}
			
//			System.out.println(multiNorms.get(columnCode));
			
			double[] vals = new double[seqs.length];
			double aa1;
			double aa2;
			
			// loop over all 3 coordinates and calculate AA distances
			for (j = 0; j < 3; j++) {
				for (int k = 0; k < seqs.length; k++){
					aa1 = rotCoords[seqs[k]][col[seqs[k]]][j];
					aa2 = rotCoords[seqs[k]][col[seqs[k]]-1][j];
					vals[k] =  coordDiffs[seqs[k]][col[seqs[k]]-1][j] = aa1-aa2;
				}
//				System.out.println(Arrays.toString(vals));
//				System.out.println("Coordinate"+(j+1)+" log likelihood:\t "+multiNorm.logDensity(vals));
				logli += multiNorm.logDensity(vals);
			}

//			for (int k = 0; k < seqs.length; k++){
//				for (j = 0; j < 3; j++) {
//					aa1 = rotCoords[seqs[k]][col[seqs[k]]][j];
//					aa2 = rotCoords[seqs[k]][col[seqs[k]]+1][j];
//					coordDiffs[seqs[k]][col[seqs[k]]][j] = aa2-aa1;
//				}
//			}
		}
		return logli;
	}
	
	
public double calcAllColumnContrib() {
		
		String[] align = curAlign;
		coordDiffs = new double[rotCoords.length][][];
		for(int i = 0; i < rotCoords.length; i++)
			coordDiffs[i] = new double[rotCoords[i].length-1][3];
		double logli = 0;
		int[] inds = new int[align.length];		// current char indices
		int[] col = new int[align.length];  
		for(int i = 0; i < align[refIndex].length(); i++) {
			for(int j = 0; j < align.length; j++)
				//Find the index of the amino acid in the sequence in the current
				//column.  Do this so you know which amino acid's coordinates to use later on
				col[j] = align[j].charAt(i) == '-' ? -1 : inds[j]++;
			double ll = columnContribBondAlign(col, align); 
			logli += ll;
		}
		return structTemp * logli;
		
	}	
	// TODO Change visibility of this to package, after moving
	// StructAlign.java to statalign.model.ext.plugins.structalign

	private boolean checkConsAlign(String[] align) {
		if(!Utils.DEBUG || curAlign == null)
			return false;
		if(align.length != curAlign.length)
			throw new Error("Inconsistency in StructAlign, alignment length: "+align.length+", "+curAlign.length);
		for(int i = 0; i < align.length; i++)
			if(!align[i].equals(curAlign[i]))
				throw new Error("Inconsistency in StructAlign, alignment: "+align[i]+", "+curAlign[i]);
		return true;
	}

	private boolean checkConsCovar(double[][] covar) {
		if(!Utils.DEBUG || fullCovar == null)
			return false;
		if(covar.length != fullCovar.length)
			throw new Error("Inconsistency in StructAlign, covar matrix length: "+covar.length+", "+fullCovar.length);
		for(int i = 0; i < covar.length; i++) {
			if(covar[i].length != fullCovar[i].length)
				throw new Error("Inconsistency in StructAlign, covar matrix "+i+" length: "+covar[i].length+", "+fullCovar[i].length);
			for(int j = 0; j < covar[i].length; j++)
				if(Math.abs(covar[i][j]-fullCovar[i][j]) > 1e-5)
					throw new Error("Inconsistency in StructAlign, covar matrix "+i+","+j+" value: "+covar[i][j]+", "+fullCovar[i][j]+", "+tau+", "+epsilon);
		}
		return true;
	}
	
	private boolean checkConsRots() {
		if(!Utils.DEBUG || rotCoords[0] == null)
			return false;
		double[][][] rots = new double[rotCoords.length][][];
		for(int i = 0; i < rots.length; i++) {
			if (rotCoords[i]==null) continue; 
			rots[i] = new double[rotCoords[i].length][];
			for(int j = 0; j < rots[i].length; j++)
				rots[i][j] = MathArrays.copyOf(rotCoords[i][j]);
		}
		calcAllRotations();
		for(int i = 0; i < rots.length; i++) {
			if (rots[i]==null) continue;
			for(int j = 0; j < rots[i].length; j++)
				for(int k = 0; k < rots[i][j].length; k++)
					if(Math.abs(rots[i][j][k]-rotCoords[i][j][k]) > 1e-5)
						throw new Error("Inconsistency in StructAlign, rotation "+i+","+j+","+k+": "+rots[i][j][k]+" vs "+rotCoords[i][j][k]);
		}
		return true;
	}

	private boolean checkConsLogLike(double logli) {
		if(!Utils.DEBUG || curLogLike == 0)
			return false;
		if(Math.abs(logli-curLogLike) > 1e-5)
			throw new Error("Inconsistency in StructAlign, log-likelihood "+logli+" vs "+curLogLike);
		return true;
	}

	
	private void addLocalEpsilonToDiagonal(double[][] subCovar, int[] notgap, int[] col) {		
		for (int i=0; i<notgap.length; i++) {
			subCovar[i][i] += Math.pow((bFactors[notgap[i]][col[notgap[i]]]),2) * epsilon / notgap.length;
		}
	}

	public class Column {
		public int[] col;
		Column(int[] x) {
			col = x.clone();
		}
		@Override
		public boolean equals(Object o) {
			if (o == this) return true;
			if (o == null || o.getClass() != this.getClass()) return false;
			
			Column x = (Column) o;
			if (x.col.length != col.length) return false;
			for (int i=0; i<x.col.length; i++) {
				if (x.col[i] != col[i]) return false;
			}
			return true;
		}
		@Override
		public int hashCode() {
			return Arrays.hashCode(col);
		}
	}
	/**
	 * extracts the specified rows and columns of a 2d array
	 * @param matrix, 2d array from which to extract; rows, rows to extract; cols, columns to extract
	 * @return submatrix
	 */
		
	private void calcAllRotations() {
		for(int i = 0; i < coords.length; i++) {
			if (coords[i]==null) continue;
			calcRotation(i);
		}		
	}
	
	public void calcRotation(int ind) {
		double[][] ci = coords[ind], rci = rotCoords[ind];
		if(rci == null)
			rci = rotCoords[ind] = new double[ci.length][];
		Rotation rot = new Rotation(new Vector3D(axes[ind]), angles[ind]);
		for(int i = 0; i < ci.length; i++) {
			rci[i] = rot.applyTo(new Vector3D(ci[i])).add(new Vector3D(xlats[ind])).toArray();
		}
	}
	// TODO Change visibility of this to package, after moving
	// StructAlign.java into statalign.model.ext.plugins.structalign.	

	
	/**
	 * return the full covariance matrix for the tree topology and branch lengths
	 */	
	public double[][] calcFullCovar(Tree tree) {
		// tree.names.length is equal to the number of vertices		
		distanceMatrix = new double[tree.names.length][tree.names.length];
		double[][] covar = new double[tree.names.length][tree.names.length];
		calcDistanceMatrix(tree.root, distanceMatrix);
		
		if (printRmsd) {
			rmsdTrace.distanceMatrix = new double[tree.names.length][tree.names.length];
			calcUnweightedDistanceMatrix(tree.root,rmsdTrace.distanceMatrix);
		}		
		
		// for hierarchical sigma, distance matrix calculation already incorporates multiplication 
		// by theta_i = sigma_i^2 / (2 tau)
		if(globalSigma){
			for(int i = 0; i < tree.names.length; i++)
				for(int j = i; j < tree.names.length; j++){
					covar[j][i] = covar[i][j] = tau * Math.exp(-linkFunction.f(distanceMatrix[i][j]) * sigma2[0] / (2*tau));
//					System.out.println(covar[j][i]);
				}
		} else{
			for(int i = 0; i < tree.names.length; i++)
				for(int j = i; j < tree.names.length; j++)
					covar[j][i] = covar[i][j] = tau * Math.exp(-distanceMatrix[i][j]);
		}
		for(int i = 0; i < tree.names.length; i++) {
			if (!localEpsilon) covar[i][i] += epsilon;			
		}
			
		
		if (localEpsilon)  multiNormsLocal = new HashMap<Column, MultiNormCholesky>();
		else multiNorms = new HashMap<Integer, MultiNormCholesky>(); 
		
		return covar;
	}	
	

	public void printTree(Vertex v, String vname){
		System.out.println(vname +"-" + v.name + ": " + v.edgeLength);
		if(v.left!=null){
			printTree(v.left, vname + "l");
			printTree(v.right, vname + "r");
		}
	}
	
	
	
	public int[] calcUnweightedDistanceMatrix(Vertex vertex, double[][] distMat){
		return calcDistanceMatrix(vertex,distMat, true);
	}
	public int[] calcDistanceMatrix(Vertex vertex, double[][] distMat){
		return calcDistanceMatrix(vertex,distMat, false);
	}
	/**
	 * recursive algorithm to traverse tree and calculate distance matrix between leaves 
	 */			
	public int[] calcDistanceMatrix(Vertex vertex, double[][] distMat, boolean unweighted){
		int[] subTree = new int[distMat.length + 1];
		
		// either both left and right are null or neither is
		if(vertex.left != null){
			int[] subLeft  = calcDistanceMatrix(vertex.left, distMat, unweighted);
			int[] subRight = calcDistanceMatrix(vertex.right, distMat, unweighted);
			int i = 0;
			while(subLeft[i] > -1){
				subTree[i] = subLeft[i];
				i++;
			}
			for(int j = 0; i+j < subTree.length; j++)
				subTree[i+j] = subRight[j];
		}
		else{
			subTree[0] = vertex.index;
			for(int j = 1; j < subTree.length; j++)
				subTree[j] = -1;
		}

		if (globalSigma || unweighted) {
			addEdgeLength(distMat, subTree, vertex.edgeLength);	
		}
		else {
			addEdgeLength(distMat, subTree, vertex.edgeLength * sigma2[vertex.index] / (2*tau));
		}
		return subTree;
	}
		
	// adds the length of the current edge to the distance between all leaves
	// of a subtree to all other leaves
	// 'rows' contains the indices of vertices in the subtree
	public void addEdgeLength(double[][] distMat, int[] subTree, double edgeLength){
		
		int i = 0;
		while(subTree[i] > -1){
			for(int j = 0; j < distMat.length; j++){  
				distMat[subTree[i]][j] += edgeLength;
				distMat[j][subTree[i]] += edgeLength;
			}
			i++;		
		}
			
		// edge length should not be added to distance between vertices in the subtree
		// subtract the value from these entries of the distance matrix
		i = 0;
		while(subTree[i] > -1){
			int j = 0;
			while(subTree[j] > -1){
				distMat[subTree[i]][subTree[j]] -= edgeLength;
				distMat[subTree[j]][subTree[i]] -= edgeLength;
				j++;
			}
			i++;
		}
	}


	@Override
	public int getParamChangeWeight() {
		// TODO test converge and tune value
		return pluginProposalWeight;
	}

	@Override
	public double logLikeModExtParamChange(Tree tree, ModelExtension ext) {
		// current log-likelihood always precomputed (regardless of whether ext == this)
		return curLogLike;
	}
	
	@Override
	public void beforeAlignChange(Tree tree, Vertex selectRoot) {
		oldAlign = curAlign;
		oldLogLi = curLogLike;
		if (Utils.DEBUG){
			System.out.println("The structure loglikelihood is:\t" + curLogLike);
		}
	}
	@Override
	public double logLikeAlignChange(Tree tree, Vertex selectRoot) {
		curAlign = tree.getState().getLeafAlign();
		curLogLike = calcAllColumnContrib();
		if (Utils.DEBUG){
			System.out.println("The structure loglikelihood is:\t" + curLogLike);
			System.out.println("coordinates:\n");
			int n = 0;
			for(double[][] k : rotCoords){
				System.out.println("Sequence "+(n++)+"\n");
				Funcs.printMatrix(k);
			}
			System.out.println("differences:\n");
			n = 0; 
			for(double[][] k : coordDiffs){
				System.out.println("Sequence "+(n++)+"\n");
				Funcs.printMatrix(k);
			}
			System.out.println("covariance matrix:\n");
			Funcs.printMatrix(fullCovar);
		}
		return curLogLike;
	}
	
	@Override
	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		curAlign = oldAlign;
		curLogLike = oldLogLi;
	}
	
	@Override
	public void beforeTreeChange(Tree tree, Vertex nephew) {
		oldDist = distanceMatrix;
		oldCovar = fullCovar;
		oldAlign = curAlign;
		oldLogLi = curLogLike;
		if (localEpsilon) oldMultiNormsLocal = multiNormsLocal;
		else oldMultiNorms = multiNorms;
	}
	@Override
	public double logLikeTreeChange(Tree tree, Vertex nephew) {		
		fullCovar = calcFullCovar(tree);
		curAlign = tree.getState().getLeafAlign();
		curLogLike = calcAllColumnContrib();
		return curLogLike;
	}
	
	@Override
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		if (showStructTree) {
			structTree.updateStructTree(tree.root,sigma2);
		}
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		distanceMatrix = oldDist;
		fullCovar = oldCovar;
		curAlign = oldAlign;
		curLogLike = oldLogLi;
		if (localEpsilon) multiNormsLocal = oldMultiNormsLocal;
		else multiNorms = oldMultiNorms;		
	}
	
	public void beforeContinuousParamChange(Tree tree) {
		//oldCovar = fullCovar;
		//oldLogLi = curLogLike;
		if (localEpsilon) oldMultiNormsLocal = multiNormsLocal;
		else oldMultiNorms = multiNorms;
	}
	public double logLikeContinuousParamChange(Tree tree) {		
		fullCovar = calcFullCovar(tree);		
		curLogLike = calcAllColumnContrib();
		return curLogLike;
	}
	public void afterContinuousParamChange(Tree tree, boolean accepted) {
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		//fullCovar = oldCovar;
		//curLogLike = oldLogLi;
		if (localEpsilon) multiNormsLocal = oldMultiNormsLocal;
		else multiNorms = oldMultiNorms;
	}
	
	@Override
	public double logLikeEdgeLenChange(Tree tree, Vertex vertex) {
		// do exactly the same as for topology change
		return logLikeTreeChange(tree, vertex);
	}
	@Override
	public void beforeEdgeLenChange(Tree tree, Vertex vertex) {
		// do exactly the same as for topology change
		beforeTreeChange(tree, vertex);
	}
	@Override
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		// do exactly the same as for topology change
		afterTreeChange(tree, vertex, accepted);
	}
	
	@Override
	public double logLikeIndelParamChange(Tree tree, Hmm hmm, McmcMove m) {
		// does not affect log-likelihood
		return curLogLike;
	}
	
	@Override
	public double logLikeSubstParamChange(Tree tree, SubstitutionModel model,
			int ind) {
		// does not affect log-likelihood
		return curLogLike;
	}
	
	@Override
	public double calcLogEm(int[] aligned) {
		return columnContrib(aligned);
	}
	
	public void connectStructTrace(StructTrace structTrace) {
		this.structTrace = structTrace;
	}
	public void connectStructTree(StructTreeVisualizer structTree) {
		this.structTree = structTree;
	}
	public void connectRmsdTrace(RmsdTrace _rmsdTrace) {
		rmsdTrace = _rmsdTrace;
		rmsdTrace.active = printRmsd;
	}
	
	/**
	 * This method is called when structures are added in the GUI.
	 */
	@Override
	public void dataAdded(File file, DataType data) {
		if(data instanceof ProteinSkeletons) {
			myButton.setEnabled(true);	
			myButton.setSelected(true);
			
			setActive(true);			
			selectable = true;
			
			// By default we'll switch on localSigma mode in GUI
			globalSigma = false;
								
			activateAssociatedPlugins();
			
			// If there is B-factor information in the .coor file
			// then we will also switch on localEpsilon by default
			if ((((ProteinSkeletons) data).bFactors) != null) {
				localEpsilon = true;
			}
		}
	}
	
	/** 
	 * Called when this plugin is activated in GUI mode, either
	 * by toggling the StructAlign button, or by reading in
	 * protein structure(s) from file.
	 */
	private void activateAssociatedPlugins() {
		structTrace.active = active;
		structTrace.postprocessWrite = active;
		structTrace.setSelected(active);
		
		// StructTree is activated by default in GUI mode, where possible
		showStructTree = !globalSigma;
		structTree.active = active & showStructTree;
		structTree.setSelected(active & showStructTree);
		structTree.postprocessWrite = active & showStructTree;
		
		// RMSDTrace is activated by default in GUI mode, but does not print to file
		rmsdTrace.active = active;
		rmsdTrace.setSelected(active);
		//rmsdTrace.postprocessWrite = active;
		rmsdTrace.postprocessWrite = false;
	}
	

	// </StructAlign>
}


