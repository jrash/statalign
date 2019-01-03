For more information, visit http://statalign.github.io
 
*** INTRODUCTION ***

StatAlign is an extendable software package for Bayesian analysis of protein, DNA and RNA sequences. Multiple alignments, phylogenetic trees and evolutionary parameters are co-estimated in a Markov Chain Monte Carlo framework, yielding posterior distributions for quantities such as tree topologies, branch lengths, alignments, and indel rates. 

Traditional methods conduct phylogenetic analysis on a single, fixed alignment, which can lead to strong biases in the inferred trees. In contrast, the joint estimation approach accounts for the interdependence between the various parameters. 

The models behind the analysis permit the comparison of evolutionarily distant sequences: the TKF92 insertion-deletion model can be coupled to a variety of different substitution models. A broad range of models for nucleotide and amino acid data is included in the package and the plug-in management system ensures that new models can be easily added.

Since joint sampling of alignments and trees is more computationally intensive than analysis with a single alignment, StatAlign does require significantly more runtime than tree-only MCMC analyses, and is best suited to the analysis of small to medium datasets (10-30 sequences). We are currently working on software developments that will extend this range of applicability.


*** USAGE ***

StatAlign is written in Java and thus requires no installation, but you must have Java 6 or newer on your system. If you do not have such a Java framework, please download the most recent one from:
  
  http://www.java.com/getjava/ 

Once you have Java, extract the StatAlign zip archive, and double-click the StatAlign.jar file. This will launch the graphical interface where sequences to analyse can be easily loaded from the menus -- see the Help page for more information.

Alternatively, StatAlign can also be run from the command line, which is recommended for automated, script-driven scenarios. To print the list of command line options, use the following command:

  java -jar StatAlign.jar -help

An example command line setup:

  java -Xmx512m -jar StatAlign.jar -ot=Fasta -mcmc=10k,100k,1k seqs.fasta

The -Xmx512m option is a standard JVM argument that sets the memory limit of the program at 512 MiBs. This is our recommended minimum, and increase it as necessary for large inputs. The -ot options selects the output alignment format, -mcmc sets MCMC parameters such as the number of burn-in steps and the number of samples, see the user documentation in the Help menu and on-line for tips on how to set these values. The input file must contain the sequences to align in Fasta format.


*** COMPILING FROM SOURCE ***

To compile (and package) StatAlign from source, you will need Gradle version 2.0 or later (see http://www.gradle.com/installation).

1. Create a runnable JAR file with the name 'StatAlign.jar':

  gradle clean build shadowJar

2. [Optional] Create a zip archive containing the newly compiled distribution (consisting of the JAR file plus examples, scripts, and documentation files), located at 'build/distributions/StatAlign-X.X.zip':

  gradle distZip


*** LICENSE ***

StatAlign is distributed under the GNU General Public License Version 3.
A copy of the license is found in LICENSE.txt.