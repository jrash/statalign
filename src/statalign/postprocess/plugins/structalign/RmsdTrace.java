package statalign.postprocess.plugins.structalign;

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import statalign.base.InputData;
import statalign.base.MCMCPars;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
import statalign.mcmc.McmcMove;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.Funcs;
import statalign.postprocess.Postprocess;
import statalign.postprocess.Track;
import statalign.postprocess.gui.StructAlignTraceGUI;
import statalign.postprocess.plugins.CurrentAlignment;
import statalign.postprocess.plugins.MpdAlignment;


public class RmsdTrace extends Postprocess {

	public int samples;
	public int cycles;

	@Override 
	public boolean createsMultipleOutputFiles() {
		return true;
	}
	
	@Override
	public ArrayList<String> getAdditionalFileExtensions() {
		ArrayList<String> result = new ArrayList<String>();
		result.add("cont.prob");
		return result;
	}
	
	public MCMCPars mcmcpars;
	
	public StructAlign structAlign;
	
	public double[][] contactProb;
	
	/** For adding to the MPD alignment panel in the GUI. */
	Track rmsdTrack = new Track(Color.RED, new double[1]);
	Track bFactorTrack = new Track(Color.GREEN, new double[1]);
	
	double maxLikelihood = Double.NEGATIVE_INFINITY;
	int sampleNumberMLE;
	String[] alignMLE, alignMLENames;
	double[] rmsdMLE, bFactorMLE;
	double epsilonMLE;
	
	public double[][] distanceMatrix;
	
	/** Determines scaling for RMSD annotation above sequence in GUI */
	double SCALE_FACTOR = 2.5;
	
	public String[] fullAlign;
	InputData input;

	public RmsdTrace() {
		screenable = false;
		outputable = true;
		postprocessable = true;
		postprocessWrite = false;
		selected = false;
		active = false;		
	}

	
	@Override
	public String getFileExtension() {
		return "rmsd";
	}

	@Override
	public void setSampling(boolean enabled) {
	}
	
	@Override
	public void init(ModelExtManager modelExtMan) {
		for(ModelExtension modExt : modelExtMan.getPluginList()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
				structAlign.connectRmsdTrace(this);
			}
		}
		// Either this is active already because we're in GUI mode, 
		// or it is not activated yet, and we only activate if printRmsd
		// is true.
		active |= structAlign.isActive() & structAlign.printRmsd;
		postprocessWrite = active & structAlign.printRmsd;
	}
	
	@Override
	public String[] getDependencies() {
		return new String[] { "statalign.postprocess.plugins.CurrentAlignment", 
							  "statalign.postprocess.plugins.MpdAlignment"};
	}

	CurrentAlignment curAli;
	MpdAlignment mpdAli;
	@Override
	public void refToDependencies(Postprocess[] plugins) {
		curAli = (CurrentAlignment) plugins[0];
		mpdAli = (MpdAlignment) plugins[1];
	}
	
	@Override
	public void beforeFirstSample(InputData inputData) {
//		for(ModelExtension modExt : getModExtPlugins()) {
//			if(modExt instanceof StructAlign) {
//				structAlign = (StructAlign) modExt;
//			}
//		}	
		contactProb = new double[structAlign.rotCoords[0].length][structAlign.rotCoords[1].length];
		
		if(!active)	return;		
		input = inputData;
		
		//double[] rad = calcGyration();
		if(postprocessWrite) {
		int leaves = structAlign.coords.length;
		try {
			outputFile.write("state\t");
/*			outputFile.write("# Each row contains pairwise mean-square deviations (msd_ij)\n");
			outputFile.write("# branch length distances (t_ij), and sequence identities (seqID_ij)\n");
			outputFile.write("# for each MCMC sample.\n");	*/		
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("msd" + i + "_" + j + "\t");
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("t" + i + "_" + j + "\t");
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("seqID" + i + "_" + j + "\t");
			
			
//			for(int i = 0; i < rad.length; i++)
//				outputFile.write(mcmc.tree.names[i] + "\t");
//			outputFile.write("\n");
//			for(int i = 0; i < rad.length; i++)
//				outputFile.write(rad[i] + "\t"); 

			outputFile.write("\n");
			outputFile.flush();
		} catch (IOException e){}
		}
		if (show) {			
			//mpdAli.addTrack(rmsdTrack);
			curAli.addTrack(rmsdTrack);
			if (structAlign.localEpsilon) curAli.addTrack(bFactorTrack);
		}			
	}
	
	@Override
	public void newPeek(State state) {
		if (!active) return;
		//if (show) {
			doUpdate(state,0);
		//}
	}
	
	private void doUpdate(State state, int sampleNumber) {				
		updateTracks(curAli.showFullAlignment ? state.getFullAlign() : state.getLeafAlign());
		if (!state.isBurnin && state.logLike > maxLikelihood) {
			maxLikelihood = state.logLike;			
			sampleNumberMLE = sampleNumber;
			alignMLE = curAli.showFullAlignment ? state.getFullAlign().clone() : state.getLeafAlign().clone();
			alignMLENames = state.name.clone();
			rmsdMLE = rmsdTrack.scores.clone();
			bFactorMLE = bFactorTrack.scores.clone();
			epsilonMLE = structAlign.epsilon;
		}
	}
	@Override
	public void newSample(State state, int no, int total) {
		samples = total;
		if(!active)
			return;
		if(postprocessWrite) {
			double[][] msd = calcMSD();
			double[][] seqID = calcSeqID();	
			calcContactProb(samples);
			try {
				outputFile.write(no*total+"\t");
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++)
						outputFile.write(msd[i][j] + "\t");
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++) 
						outputFile.write(distanceMatrix[i][j] + "\t");
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++)
						outputFile.write(seqID[i][j] + "\t");
				outputFile.write("\n");
			} catch (IOException e){
				e.printStackTrace();
			}		
		}		
		//if (show) {
			doUpdate(state,no);
		//}
	}
	
	@Override
	public void afterLastSample() {
		if(!active)
			return;		
		if(postprocessWrite) {
		try {
			if(structAlign.rotCoords.length == 2){
				for(int i=0; i<contactProb.length; i++){
					for(int j=0; j<contactProb[0].length; j++)
						contactProb[i][j] /= samples;
				}
				printMatrix(contactProb, additionalOutputFiles.get(0));
			}
			outputFile.close();
			additionalOutputFiles.get(0).close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}				
		try { 
			if (rmsdMLE != null) {
				FileWriter mle = new FileWriter(getBaseFileName()+"mle." + getFileExtension());
				mle.write("# Maximum likelihood = "+maxLikelihood+" at sample "+sampleNumberMLE+"\n");
				mle.write("# RMSD\tAverage B-factor\n");
				mle.flush();
				for (int i=0; i<rmsdMLE.length; i++) {
					boolean allGap = true;					
					for (int j=0; j<structAlign.rotCoords.length; j++) { // Assume first sequences are non-internals
						if (alignMLE[j].charAt(i) != '-') allGap = false;
					}
					if (!allGap) {
						mle.write(rmsdMLE[i]+"");
						if (structAlign.localEpsilon) mle.write("\t"+3*epsilonMLE*bFactorMLE[i]);
						mle.write("\n");
					}
					
				}
				mle.close();
				
				FileWriter aliMLE = new FileWriter(getBaseFileName()+"mle.ali");

				// Form an array from the leaves of the alignment 
				String[] aln = Utils.alignmentTransformation(alignMLE, alignMLENames,
						"Fasta", input);
				for (int i = 0; i < aln.length; i++) {
					aliMLE.write(aln[i] + "\n");
				}
				aliMLE.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		}
	}
	
	public static void printMatrix(double[][] m) {
		for(int i = 0; i < m.length; i++)
			System.out.println(Arrays.toString(m[i]));
		System.out.println();
	}
	
	public static void printMatrix(double[][] m, FileWriter outputFile) throws IOException {
		for(int i = 0; i < m.length; i++)
			outputFile.write(Arrays.toString(m[i])+"\n");
		outputFile.write("\n");
	}
	
	public double[][] calcMSD(){
		double[][][] coor = structAlign.rotCoords;
		String[] align = structAlign.curAlign;
		if (coor.length == 2)
			for(int i=0; i< coor.length;i++){
				RealMatrix temp = new Array2DRowRealMatrix(coor[i]);
				RealVector mean = Funcs.meanVectorAligned(temp,align,i);
				for(int j = 0; j < coor[i].length; j++)
					 coor[i][j]= temp.getRowVector(j).subtract(mean).toArray();
			}
		else
			for(int i=0; i< coor.length;i++){
				RealMatrix temp = new Array2DRowRealMatrix(coor[i]);
				RealVector mean = Funcs.meanVector(temp);
				for(int j = 0; j < coor[i].length; j++)
					coor[i][j]= temp.getRowVector(j).subtract(mean).toArray();
			}
		int leaves = coor.length;
		boolean[] hasStructure = new boolean[leaves];		
		boolean igap, jgap;
		double[][] msd = new double[leaves][leaves];
		for (int i=0; i<leaves; i++) {
			hasStructure[i] = (structAlign.coords[i]!=null);				
		}
		for(int i = 0; i < leaves-1; i++){		
			for(int j = i+1; j < leaves; j++){
				int ii = 0, jj = 0, n = 0;
				for(int k = 0; k < align[0].length(); k++){
					igap = align[i].charAt(k) == '-';
					jgap = align[j].charAt(k) == '-';
					if(!igap & !jgap & hasStructure[i] & hasStructure[j]){
						msd[i][j] += sqDistance(coor[i][ii], coor[j][jj]);
						n++;
					}
					ii += igap ? 0 : 1;
					jj += jgap ? 0 : 1;
				}
				msd[i][j] /= n;
			}
		}
		return msd;	}
	public void updateTracks(String[] align){
		double[][][] coor = structAlign.rotCoords;
		int leaves = coor.length;
		boolean igap, jgap;		
		//int n = leaves * (leaves-1) / 2; // Number of pairwise comparisons
		
		int alignmentLength = align[0].length();
		int leafAlignmentLength = align[0].length();
		boolean[] allGapped = new boolean[alignmentLength];
		boolean[] allGappedLeaf = new boolean[alignmentLength];
		for(int k = 0; k < align[0].length(); k++){
			allGapped[k] = true;
			for(int i = 0; i < align.length; i++){
				allGapped[k] &= (align[i].charAt(k) == '-');				
				if (i < leaves) allGappedLeaf[k] &= (align[i].charAt(k) == '-');
			}
			if (allGapped[k]) alignmentLength--;
			if (allGappedLeaf[k]) leafAlignmentLength--;
		}
		rmsdTrack.scores = new double[alignmentLength];
		rmsdTrack.max = 0.0; rmsdTrack.min = Double.POSITIVE_INFINITY; rmsdTrack.mean = 0.0;
		if (structAlign.localEpsilon) {
			bFactorTrack.scores = new double[alignmentLength];
			bFactorTrack.max = 0.0; bFactorTrack.min = Double.POSITIVE_INFINITY; bFactorTrack.mean = 0.0;
		}
		int[] index = new int[leaves];		
		for(int k = 0, kk=0; k < align[0].length(); k++){				

			int n=0;
			if (allGapped[k]) continue;			
			for(int i = 0; i < leaves-1; i++){
				for(int j = i+1; j < leaves; j++){													
					igap = align[i].charAt(k) == '-';
					jgap = align[j].charAt(k) == '-';
					if(!igap & !jgap & !(coor[i]==null) & !(coor[j]==null)){
						rmsdTrack.scores[kk] += sqDistance(coor[i][index[i]], coor[j][index[j]]);
						++n;
					}					
				}				
			}			
			if (n > 0) rmsdTrack.scores[kk] /= n;
			if (n==0)  rmsdTrack.scores[kk] = Double.NaN;
			if (!Double.isNaN(rmsdTrack.scores[kk])) {			
				if (rmsdTrack.scores[kk] > 0) rmsdTrack.scores[kk] = Math.sqrt(rmsdTrack.scores[kk]);
				
				if (rmsdTrack.scores[kk] > rmsdTrack.max) rmsdTrack.max = rmsdTrack.scores[kk];
				if (rmsdTrack.scores[kk] > 0 && rmsdTrack.scores[kk] < rmsdTrack.min) rmsdTrack.min = rmsdTrack.scores[kk];
				rmsdTrack.mean += rmsdTrack.scores[kk] / leafAlignmentLength;
			}
			
			int nBfactor = 0;
			for(int i = 0; i < leaves; i++) {
				if (align[i].charAt(k) != '-') {
					if(structAlign.coords[i]==null) continue;
					if (structAlign.localEpsilon) bFactorTrack.scores[kk] += Math.pow(structAlign.bFactors[i][index[i]],2) * structAlign.epsilon;
					nBfactor++;
					index[i]++;
				}
			}			 	
			if (structAlign.localEpsilon) {
				bFactorTrack.scores[kk] /= nBfactor;
				bFactorTrack.scores[kk] = Math.sqrt(bFactorTrack.scores[kk]);				 
				if (bFactorTrack.scores[kk] == 0) bFactorTrack.scores[kk] = Double.NaN; 									
				if (!Double.isNaN(bFactorTrack.scores[kk])) { 							
					if (bFactorTrack.scores[kk] > bFactorTrack.max) bFactorTrack.max = bFactorTrack.scores[kk];
					if (bFactorTrack.scores[kk] > 0 && bFactorTrack.scores[kk] < bFactorTrack.min) bFactorTrack.min = bFactorTrack.scores[kk];
					bFactorTrack.mean += bFactorTrack.scores[kk] / leafAlignmentLength;		
				}					
			}
			
			kk++;			
		}				
	}
	
	public double sqDistance(double[] x, double[] y){		
		double d = 0;
		for(int i = 0; i < x.length; i++)
			d += Math.pow(x[i] - y[i], 2.0);
		return d;
	}
	
	public double[] calcGyration(){
		double[][][] coor = structAlign.coords;
		int leaves = coor.length;
		double[] radii = new double[leaves];
		for(int i = 0; i < leaves; i++){
			radii[i] = 0;
			// coordinates are centered in StructAlign.initRun()
			for(int j = 0; j < coor[i].length; j++)
				for(int k = 0; k < coor[i][0].length; k++)
					radii[i] += Math.pow(coor[i][j][k], 2.0);
			radii[i] /= coor[0].length;
			radii[i] = Math.pow(radii[i], 0.5);
		}
		return radii;
	}
	
	public double[][] calcSeqID(){
		String[] align = structAlign.curAlign;
		int leaves = align.length;
		boolean igap, jgap;
		double[][] seqID = new double[leaves][leaves];
		for(int i = 0; i < leaves-1; i++){
			for(int j = i+1; j < leaves; j++){
				double match = 0, id = 0;
				for(int k = 0; k < align[0].length(); k++){
					igap = align[i].charAt(k) == '-';
					jgap = align[j].charAt(k) == '-';
					if(!igap & !jgap){
						id += align[i].charAt(k) == align[j].charAt(k) ? 1 : 0;
						match++;
					}
				}
				seqID[i][j] = id / match;
			}
		}
		return seqID;
	}
	
	public void calcContactProb(int total){
		String[] align = structAlign.curAlign;
		boolean ngap, mgap;
		int aaIdxN = 0, aaIdxM = 0;
		for(int k = 0; k < align[0].length(); k++){
			ngap = align[0].charAt(k) == '-';
			mgap = align[1].charAt(k) == '-';
			if(!ngap & !mgap){
				contactProb[aaIdxN][aaIdxM]++;
			}
			if(!ngap)
				aaIdxN++;
			if(!mgap)
				aaIdxM++;
		}
	}


	@Override
	public String getTabName() {
		// TODO Auto-generated method stub
		return "RMSD trace";
	}


	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public JPanel getJPanel() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public String getTip() {
		// TODO Auto-generated method stub
		return "";
	}
}
