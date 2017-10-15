import NumericalMethod.HLL;
import NumericalMethod.LxF;
import java.io.*;

/*
	This java program is a final year project in 
	CUHK Physics surpervised by Dr Leung Po Kin.
	The codes were developed under mainly imperative
	ways. Although I have tried to implement in OOP 
	styles, the output of the program still needs the
	users to implements. 
	
	This project mainly focused on HLL method. LxF method
	had also implemented but not well-organized.

	HLL has already passed the shock tube test. Due to 
	the time limitation, Kelvin-Helmholtz instability test
	had not been studies well.

	Note: LxF class is not well-organized. It will be hard 
	to read and used. Please focus on HLL class if you are
	interested.
*/

class Main {
	public static void main(String[] args) throws IOException, InterruptedException {
		// Provide the number of grids, size of box and target test to constructor
		HLL HLL = new HLL(250,1,"1D_y_shock_tube");
		// This method will write the output into files
		HLL.iterate(0.2);
	}
}