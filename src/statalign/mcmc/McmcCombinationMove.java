package statalign.mcmc;

import java.util.List;
import java.util.ArrayList;
import static statalign.utils.Reversed.reversed;

public class McmcCombinationMove extends McmcMove {

	protected List<McmcMove> mcmcMoves = new ArrayList<McmcMove>();
	
	public McmcCombinationMove (ArrayList<McmcMove> moves) {
		if (moves.size() < 2) {
			throw new IllegalArgumentException("McmcCombinationMove must contain at least two McmcMove objects");
		}		
		name = moves.get(0).name;
		mcmcMoves.add(moves.get(0));
		for (int i=1; i<moves.size(); i++) {
			name += "_"+moves.get(i).name;
			mcmcMoves.add(moves.get(i));
		}
	}
	public void copyState(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.copyState(externalState);
		}
	}
	public double proposal(Object externalState) {
		double logProposalRatio = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			double oldProposalWidthControlVariable = mcmcMove.proposalWidthControlVariable;
			mcmcMove.proposalWidthControlVariable *= proposalWidthControlVariable;
			logProposalRatio += mcmcMove.proposal(externalState); 
			mcmcMove.proposalWidthControlVariable = oldProposalWidthControlVariable;
		}
		// We implicitly assume here that the order of the proposals
		// does not matter.
		return logProposalRatio;
	}
	
	@Override
	public McmcModule getOwner() { 
		return mcmcMoves.get(0).getOwner();
		// NB it doesn't matter which McmcMove we use here, since it is
		// only used to call back to the Mcmc object running
		// all of them, but since the McmcCombinationMove does not necessarily
		// have a unique ModelExtension as its owner, we must use one of
		// its McmcMove objects to do the callback.
	}
			
	public double logPriorDensity(Object externalState) {
		double logPriorDensity = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			logPriorDensity += mcmcMove.logPriorDensity(externalState);
		}
		return logPriorDensity;
	}
	public void updateLikelihood(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.updateLikelihood(externalState);
		}
	}
	public void afterMove(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.afterMove(externalState);
		}
	}
	public void restoreState(Object externalState) {
		for (McmcMove mcmcMove : reversed(mcmcMoves)) {
			mcmcMove.restoreState(externalState);
		}
	}
}
