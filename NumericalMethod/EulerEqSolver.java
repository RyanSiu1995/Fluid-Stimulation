package NumericalMethod;

interface EulerEqSolver {
	/* 
		Required Qualites:
		Gamma = Heat Capacity Ratio
		CFL = CFL condition
		L is the length of the box in arbitary unit
	*/
    public double gamma = 1.4;
    public double CFL = 0.7;
    public double L = 1;
    /*
		Required Method:
		setIC = set the initial conditions
		updateFlux = update the flux in 
		the edges of boxes
		updateState = update all the qualities 
		in next time step
		checkCFL = check the CFL condition to
		change the time step
    */
    public void setIC();
    public void updateFlux();
    public void updateState();
    public void checkCFL();
}