package statalign.model.ext.plugins.structalign;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.model.ext.plugins.structalign.Transformation;

public class Funcs{
	
	public static double[][] getSubMatrix(double[][] matrix, int[] rows, int[] cols) {
		double[][] submat = new double[rows.length][cols.length];
		for(int i = 0; i < rows.length; i++)
			for(int j = 0; j < cols.length; j++)
				submat[i][j] = matrix[rows[i]][cols[j]];
		return submat;
	}
	/** Calculates the cross product of 2 vectors
	 * 
	 * @param a first vector
	 * @param b second vector
	 * @return cross product
	 */
	/*public static RealVector crossProduct(RealVector a, RealVector b){
		RealMatrix skew = new Array2DRowRealMatrix(
				new double[][] {{0, -a.getEntry(2), a.getEntry(1)}, {a.getEntry(2), 0, -a.getEntry(0)},
						{-a.getEntry(1), a.getEntry(0), 0}});
		return skew.operate(b);
	}		*/

	/** 
	 * Calculates the rotation matrix from an axis and angle
	 * with axis x and angle r, the rotation matrix is
	 * cos(r)I + sin(r)cross(x) + (1-cos(r))outer(x)
	 * where cross(x) is the crossproduct matrix of x and outer(x) is the outer (tensor) product
	 * 
	 * @return Transformation object with rotation matrix calculated from axis and angle
	 */
	/*public static RealMatrix calcRotationMatrix(RealVector axis, double rot){
		RealMatrix outer = axis.outerProduct(axis);
		// use transpose of cross product matrix (as given on wikipedia) because
		// we post-multiply by rotation matrix (other components of final step are symmetric)
		RealMatrix crossTranspose = new Array2DRowRealMatrix(new double[][] { 
				{0, axis.getEntry(2), -axis.getEntry(1)},
				{-axis.getEntry(2), 0, axis.getEntry(0)},
				{axis.getEntry(1), -axis.getEntry(0), 0} 
		});
		RealMatrix Icos = new Array2DRowRealMatrix(new double[3][3]);
		for(int i = 0; i < 3; i++)
			Icos.setEntry(i, i, Math.cos(rot));
		crossTranspose = crossTranspose.scalarMultiply(Math.sin(rot));
		outer = outer.scalarMultiply(1 - Math.cos(rot));
		return Icos.add(crossTranspose).add(outer);
	} */

	/**
	 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector
	 * @param A - coordinate matrix
	 * @return mean vector
	 */
	
	public static RealVector meanVector(RealMatrix A){
		RealVector mean = new ArrayRealVector(new double[3]);
		for(int i = 0; i < 3; i ++){
			for(int j = 0; j < A.getColumn(0).length; j++)
				mean.addToEntry(i, A.getEntry(j, i));
			mean.setEntry(i, mean.getEntry(i) / A.getColumn(0).length);
		}
		return mean;
	}


	public static Vertex sampleVertex(Tree tree){
		/* Vertex v;
		int n = tree.vertex.length;		// number of vertices
		int l = coords.length;			// number of leaves
		
		
		List<Integer> inds = findRefSubtrees(tree, 0);	// returns indices of all ancestor vertices of reference protein
		if(inds.size() + l < n){
			int prop = Utils.generator.nextInt(n-l) + l;		// don't choose a leaf vertex			
			while(inds.contains(new Integer(prop)))
				prop = Utils.generator.nextInt(n-l) + l;		// don't choose a subtree containing reference protein
			v = tree.vertex[prop];
			return v;
		} else
			return tree.root; */
		Vertex v = tree.vertex[Utils.generator.nextInt(tree.vertex.length)];
		while (v == tree.root) {
			v = tree.vertex[Utils.generator.nextInt(tree.vertex.length)];
		}
		return  v;
	}
	
	public static ArrayList<Integer> findRefSubtrees(Tree tree, int refInd){
		ArrayList<Integer> inds = new ArrayList<Integer>(0);
		moveUp(tree.vertex[refInd].parent, inds);
		return inds;
	}
	
	public static void moveUp(Vertex v, List<Integer> inds){
		inds.add(v.index);
		if(v.parent != null)
			moveUp(v.parent, inds);
	}
	
	public static ArrayList<Integer> collectLeaves(Vertex v){
		ArrayList<Integer> inds = new ArrayList<Integer>(0);
		moveDown(v, inds);
		return inds;
	}
	
	public static void moveDown(Vertex v, List<Integer> inds){
		if(v.left != null){
			moveDown(v.left, inds);
			moveDown(v.right, inds);
		}
		else
			inds.add(v.index);
	}
	
	public static void printMatrix(double[][] ls) {

		for (int i = 0; i < ls.length; i++) {
			for (int j = 0; j < ls[i].length; j++) {
				System.out.print(ls[i][j] + "\t");
			}
			System.out.println();
		}

	}
	
	public static void writeRotationFiles(Tree tree, double[][][] coords, 
			double[][] xlats, double[][] axes, double[] angles){
		try{
			FileWriter fstream = new FileWriter("names.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for(int i = 0; i < coords.length; i++)
				out.write(tree.vertex[i].name + "\n");
			out.close();
	
			FileWriter fstream2 = new FileWriter("trans.txt");
			BufferedWriter out2 = new BufferedWriter(fstream2);
			for(int i = 0; i < coords.length; i++){
				for(int j = 0; j < 3; j++)
					out2.write(xlats[i][j] + "\t");
				out2.write("\n");
			}
			out2.close();
			
			FileWriter fstream3 = new FileWriter("rots.txt");
			BufferedWriter out3 = new BufferedWriter(fstream3);
			for(int i = 0; i < coords.length; i++){
				if (coords[i]==null) continue;
				Transformation trans = new Transformation(axes[i], angles[i], xlats[i]);
				RealMatrix temp = trans.getRealRotation();
				for(int j = 0; j < 3; j++){
					for(int k = 0; k < 3; k++)
						out3.write(temp.getEntry(j, k) + "\t");
					out3.write("\n");
				}
			}
			out3.close();
			
		} catch (Exception e){
			System.err.println("File writing error: " + e.getMessage());
		}
		
	}
	
	public static void initLSRotations(Tree tree, double[][][] coords, 
			double[][] xlats, double[][] axes, double[] angles){
		String[] align = tree.getState().getLeafAlign();		
		int refIndex = 0;
		while (coords[refIndex]==null) ++refIndex;		
		String ref = align[refIndex];
		axes[refIndex] = new double[] { 1, 0, 0 };
		angles[refIndex] = 0;
		xlats[refIndex] = new double[] { 0, 0, 0 };
		for(int i = 0; i < align.length; i++){
			if (i==refIndex || coords[i]==null) continue; 
			ArrayList<Integer> refInds = new ArrayList<Integer>(0);
			ArrayList<Integer> otherInds = new ArrayList<Integer>(0);
			String other = align[i];
			
			int r = 0, o = 0;
			for(int j = 0; j < align[refIndex].length(); j++){
				if(ref.charAt(j) != '-' & other.charAt(j) != '-'){
					refInds.add(r);
					otherInds.add(o);
				}
				r += ref.charAt(j) != '-' ? 1 : 0;
				o += other.charAt(j) != '-' ? 1 : 0;
			}
			
			double[][] refSub = getRowSub(coords[refIndex], refInds);
			double[][] otherSub = getRowSub(coords[i], otherInds);
			
			Transformation trans = new Transformation();
			
			trans.lsTrans(new Array2DRowRealMatrix(refSub), new Array2DRowRealMatrix(otherSub));
			axes[i] = trans.rotMatrix.getAxis().toArray();
			xlats[i] = trans.xlat.toArray();
			angles[i]  = trans.rotMatrix.getAngle();
		}
	}
	
	public static double[][] getRowSub(double[][] coord, ArrayList<Integer> rows){
		double[][] sub = new double[rows.size()][coord[0].length];
		
		for(int i = 0; i < sub.length; i++)
			for(int j = 0; j < sub[0].length; j++)
				sub[i][j] = coord[rows.get(i)][j];
		return sub;
	}

	/**
	 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector using only aligned positions
	 * @param A - coordinate matrix
	 * @return mean vector
	 */
	
	//Test this
	
	public static RealVector meanVectorAligned(RealMatrix A, String[] align, int seq){
		
		RealVector mean = new ArrayRealVector(new double[3]);
		
		boolean ngap, mgap, curSeqGap;
		int aaIdx = 0, aligned = 0;
		for(int j = 0; j < align[0].length(); j++){
			ngap = align[0].charAt(j) == '-';
			mgap = align[1].charAt(j) == '-';
			curSeqGap = align[seq].charAt(j) == '-';
			if(!ngap & !mgap){
				aligned++;
				for(int i = 0; i < 3; i ++){
					mean.addToEntry(i, A.getEntry(aaIdx,i));
				}
			}	
			if(!curSeqGap)
				aaIdx++;
		}
		
		for(int i = 0; i < 3; i ++){
			mean.setEntry(i, mean.getEntry(i) / aligned);
		}
		return mean;
	}

	
	
	
}


