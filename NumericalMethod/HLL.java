package NumericalMethod;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class HLL implements EulerEqSolver {
	//state, pressure can be calculated directly
    private double rho[][];
    private double rho_u[][];
    private double rho_v[][];
    private double rho_w[][];
    private double pressure[][];
    private double E[][];
    //flux, no pressure flux exist! Note that (0,0) is useless
    private double x_rhoF[][];
    private double x_rho_uF[][];
    private double x_rho_vF[][];
    private double x_rho_wF[][];
    private double x_EF[][];
    private double y_rhoF[][];
    private double y_rho_uF[][];
    private double y_rho_vF[][];
    private double y_rho_wF[][];
    private double y_EF[][];
    //Cell related parameter
    private int size, ghost;
    private double dx, dy, dt;
    private String name;

    //Inititalize the number of cell
    public HLL(int size, int ghostcell, String name) {
        this.size = size;
        this.ghost = ghostcell;
        this.name = name;
        dx = L/(double)size;
        dy = dx;
        rho = new double[size+ghostcell*2][size+ghostcell*2];
        rho_u = new double[size+ghostcell*2][size+ghostcell*2];
        rho_v = new double[size+ghostcell*2][size+ghostcell*2];
        rho_w = new double[size+ghostcell*2][size+ghostcell*2];
        pressure = new double[size+ghostcell*2][size+ghostcell*2];
        E = new double[size+ghostcell*2][size+ghostcell*2];
        x_rhoF = new double[size+ghostcell][size+ghostcell];
        x_rho_uF= new double[size+ghostcell][size+ghostcell];
        x_rho_vF= new double[size+ghostcell][size+ghostcell];
        x_rho_wF= new double[size+ghostcell][size+ghostcell];
        x_EF = new double[size+ghostcell][size+ghostcell];
        y_rhoF= new double[size+ghostcell][size+ghostcell];
        y_rho_uF= new double[size+ghostcell][size+ghostcell];
        y_rho_vF= new double[size+ghostcell][size+ghostcell];
        y_rho_wF= new double[size+ghostcell][size+ghostcell];
        y_EF = new double[size+ghostcell][size+ghostcell];
    }

    //Initialize the initial condition
    @Override
    public void setIC() {
    	// In order to carry out the testes
    	// Switch is used here
        switch(name) {
        	// 2D shock tube
            case "2D shock tube":
                for(int i = ghost; i < size+ghost; i++) {
                    for(int j = ghost; j < size+ghost; j++) {
                        if(i > j) {
                            rho[i][j] = 1;
                            pressure[i][j] = 1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else {
                            rho[i][j] = 0.125;
                            pressure[i][j] = 0.1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                    }
                }
                setBoundary(1);
                checkCFL();
                System.out.println("Finish setIC");
                break;
            case "1D_x_shock_tube":
                // 1D shock tube
                for(int i = ghost; i < size+ghost; i++) {
                    for(int j = ghost; j < size+ghost; j++) {
                        if(i < size/2) {
                            rho[i][j] = 1;
                            pressure[i][j] = 1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else {
                            rho[i][j] = 0.125;
                            pressure[i][j] = 0.1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                    }
                }
                setBoundary(2);
                checkCFL();
                System.out.println("Finish setIC");
                break;
            case "1D_y_shock_tube":
                // 1D shock tube
                for(int i = ghost; i < size+ghost; i++) {
                    for(int j = ghost; j < size+ghost; j++) {
                        if(j < size / 2){
                            rho[i][j] = 1;
                            pressure[i][j] = 1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else {
                            rho[i][j] = 0.125;
                            pressure[i][j] = 0.1;
                            rho_u[i][j] = 0;
                            rho_v[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                    }
                }
                setBoundary(3);
                checkCFL();
                System.out.println("Finish setIC");
                break;
            // This test cannot be passed
            case "Kelvin-Helmholtz Instability Test":
                for(int i = ghost; i < size+ghost; i++) {
                    for(int j = ghost; j < size+ghost; j++) {
                        if(i < size/3) {
                            rho[i][j] = 1;
                            pressure[i][j] = 2.5;
                            rho_v[i][j] = -0.5*1;
                            rho_u[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else 
                        if(i == size/3) {
                            rho[i][j] = 2;
                            pressure[i][j] = 2.5;
                            rho_v[i][j] = 0.5*2;
                            rho_u[i][j] = rho[i][j] * Math.sin(2*Math.PI*(double) j/ ((double)size/2))/100;
                            //rho_u[i][j] = 2*Math.sin(2 * Math.PI * (double)j / (double)32)/100;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else
                        if(i < 2*size/3) {
                            rho[i][j] = 2;
                            pressure[i][j] = 2.5;
                            rho_v[i][j] = 0.5*2;
                            rho_u[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else
                        if(i == 2*size/3) {
                            rho[i][j] = 2;
                            pressure[i][j] = 2.5;
                            rho_v[i][j] = 0.5*2;
                            rho_u[i][j] = rho[i][j] * Math.sin(2*Math.PI*(double) j/ ((double)size/2))/100;
                            //rho_u[i][j] = 2*Math.sin((double)i/(double)32 * 2 *Math.PI)/100;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                        else {
                            rho[i][j] = 1;
                            pressure[i][j] = 2.5;
                            rho_v[i][j] = -0.5*1;
                            rho_u[i][j] = 0;
                            rho_w[i][j] = 0;
                            E[i][j] = E_finding(i, j);
                        }
                    }
                }
                setBoundary(4);
                checkCFL();
                System.out.println("Finish setIC");
                break;
        }
    }

    // Set the Boundarg Conditions for different cases
    private void setBoundary(int a){
        switch(a){
            case 1: 
            	// All absorpt boundary
                for(int i = ghost;i <= size+ghost; i++) {
                    rho[i][0] = rho[i][1];
                    pressure[i][0] = pressure[i][1];
                    rho_u[i][0] = rho_u[i][1];
                    rho_v[i][0] = rho_v[i][1];
                    rho_w[i][0] = rho_w[i][1];
                    E[i][0] = E[i][1];
                    rho[i][size+ghost] = rho[i][size-1+ghost];
                    pressure[i][size+ghost] = pressure[i][size-1+ghost];
                    rho_u[i][size+ghost] = rho_u[i][size-1+ghost];
                    rho_v[i][size+ghost] = rho_v[i][size-1+ghost];
                    rho_w[i][size+ghost] = rho_w[i][size-1+ghost];
                    E[i][size+ghost] = E[i][size-1+ghost];
                }
                for(int j = ghost;j <= size+ghost; j++) {
                    rho[0][j] = rho[1][j];
                    pressure[0][j] = pressure[1][j];
                    rho_u[0][j] = rho_u[1][j];
                    rho_v[0][j] = rho_v[1][j];
                    rho_w[0][j] = rho_w[1][j];
                    E[0][j] = E[1][j];
                    rho[size+ghost][j] = rho[size+ghost-1][j];
                    pressure[size+ghost][j] = pressure[size+ghost-1][j];
                    rho_u[size+ghost][j] = rho_u[size+ghost-1][j];
                    rho_v[size+ghost][j] = rho_v[size+ghost-1][j];
                    rho_w[size+ghost][j] = rho_w[size+ghost-1][j];
                    E[size+ghost][j] = E[size+ghost-1][j];
                }
                //Fix four conner 
                rho[0][0] = rho[1][0];
                rho[rho.length-1][0] = rho[rho.length-1][1];
                rho[0][rho.length-1] = rho[1][rho.length-1];
                rho[rho.length-1][rho.length-1] = rho[rho.length-2][rho.length-1];
                E[0][0] = E[1][0];
                E[rho.length-1][0] = E[rho.length-1][1];
                E[0][rho.length-1] = E[1][rho.length-1];
                E[rho.length-1][rho.length-1] = E[rho.length-2][rho.length-1];
                rho_u[0][0] = rho_u[1][0];
                rho_u[rho.length-1][0] = rho_u[rho.length-1][1];
                rho_u[0][rho.length-1] = rho_u[1][rho.length-1];
                rho_u[rho.length-1][rho.length-1] = rho_u[rho.length-2][rho.length-1];
                rho_v[0][0] = rho_v[1][0];
                rho_v[rho.length-1][0] = rho_v[rho.length-1][1];
                rho_v[0][rho.length-1] = rho_v[1][rho.length-1];
                rho_v[rho.length-1][rho.length-1] = rho_v[rho.length-2][rho.length-1];
                rho_w[0][0] = rho_w[1][0];
                rho_w[rho.length-1][0] = rho_w[rho.length-1][1];
                rho_w[0][rho.length-1] = rho_w[1][rho.length-1];
                rho_w[rho.length-1][rho.length-1] = rho_w[rho.length-2][rho.length-1];
                pressure[0][0] = pressure[1][0];
                pressure[rho.length-1][0] = pressure[rho.length-1][1];
                pressure[0][rho.length-1] = pressure[1][rho.length-1];
                pressure[rho.length-1][rho.length-1] = pressure[rho.length-2][rho.length-1];
                break;
            case 2:
            	// y out flow, x periodic
                for(int i = ghost;i <= size+ghost; i++) {
                    rho[i][0] = rho[i][size+ghost-1];
                    pressure[i][0] = pressure[i][size+ghost-1];
                    rho_u[i][0] = rho_u[i][size+ghost-1];
                    rho_v[i][0] = rho_v[i][size+ghost-1];
                    rho_w[i][0] = rho_w[i][size+ghost-1];
                    E[i][0] = E[i][size+ghost-1];
                    rho[i][size+ghost] = rho[i][ghost];
                    pressure[i][size+ghost] = pressure[i][ghost];
                    rho_u[i][size+ghost] = rho_u[i][ghost];
                    rho_v[i][size+ghost] = rho_v[i][ghost];
                    rho_w[i][size+ghost] = rho_w[i][ghost];
                    E[i][size+ghost] = E[i][ghost];
                }
                for(int j = ghost;j <= size+ghost; j++) {
                    rho[0][j] = rho[1][j];
                    pressure[0][j] = pressure[1][j];
                    rho_u[0][j] = rho_u[1][j];
                    rho_v[0][j] = rho_v[1][j];
                    rho_w[0][j] = rho_w[1][j];
                    E[0][j] = E[1][j];
                    rho[size+ghost][j] = rho[size+ghost-1][j];
                    pressure[size+ghost][j] = pressure[size+ghost-1][j];
                    rho_u[size+ghost][j] = rho_u[size+ghost-1][j];
                    rho_v[size+ghost][j] = rho_v[size+ghost-1][j];
                    rho_w[size+ghost][j] = rho_w[size+ghost-1][j];
                    E[size+ghost][j] = E[size+ghost-1][j];
                }
                //Fix four conner 
                rho[0][0] = rho[1][0];
                rho[rho.length-1][0] = rho[rho.length-1][1];
                rho[0][rho.length-1] = rho[1][rho.length-1];
                rho[rho.length-1][rho.length-1] = rho[rho.length-2][rho.length-1];
                E[0][0] = E[1][0];
                E[rho.length-1][0] = E[rho.length-1][1];
                E[0][rho.length-1] = E[1][rho.length-1];
                E[rho.length-1][rho.length-1] = E[rho.length-2][rho.length-1];
                rho_u[0][0] = rho_u[1][0];
                rho_u[rho.length-1][0] = rho_u[rho.length-1][1];
                rho_u[0][rho.length-1] = rho_u[1][rho.length-1];
                rho_u[rho.length-1][rho.length-1] = rho_u[rho.length-2][rho.length-1];
                rho_v[0][0] = rho_v[1][0];
                rho_v[rho.length-1][0] = rho_v[rho.length-1][1];
                rho_v[0][rho.length-1] = rho_v[1][rho.length-1];
                rho_v[rho.length-1][rho.length-1] = rho_v[rho.length-2][rho.length-1];
                rho_w[0][0] = rho_w[1][0];
                rho_w[rho.length-1][0] = rho_w[rho.length-1][1];
                rho_w[0][rho.length-1] = rho_w[1][rho.length-1];
                rho_w[rho.length-1][rho.length-1] = rho_w[rho.length-2][rho.length-1];
                pressure[0][0] = pressure[1][0];
                pressure[rho.length-1][0] = pressure[rho.length-1][1];
                pressure[0][rho.length-1] = pressure[1][rho.length-1];
                pressure[rho.length-1][rho.length-1] = pressure[rho.length-2][rho.length-1];
                break;
            case 3: 
            	//x out flow, y periodic
                for(int i = ghost;i <= size+ghost; i++) {
                    rho[0][i] = rho[size+ghost-1][i];
                    pressure[0][i] = pressure[size+ghost-1][i];
                    rho_u[0][i] = rho_u[size+ghost-1][i];
                    rho_v[0][i] = rho_v[size+ghost-1][i];
                    rho_w[0][i] = rho_w[size+ghost-1][i];
                    E[0][i] = E[size+ghost-1][i];
                    rho[size+ghost][i] = rho[ghost][i];
                    pressure[size+ghost][i] = pressure[ghost][i];
                    rho_u[size+ghost][i] = rho_u[ghost][i];
                    rho_v[size+ghost][i] = rho_v[ghost][i];
                    rho_w[size+ghost][i] = rho_w[ghost][i];
                    E[size+ghost][i] = E[ghost][i];
                }
                for(int j = ghost;j <= size+ghost; j++) {
                    rho[j][0] = rho[j][1];
                    pressure[j][0] = pressure[j][1];
                    rho_u[j][0] = rho_u[j][1];
                    rho_v[j][0] = rho_v[j][1];
                    rho_w[j][0] = rho_w[j][1];
                    E[j][0] = E[j][1];
                    rho[j][size+ghost] = rho[j][size+ghost-1];
                    pressure[j][size+ghost] = pressure[j][size+ghost-1];
                    rho_u[j][size+ghost] = rho_u[j][size+ghost-1];
                    rho_v[j][size+ghost] = rho_v[j][size+ghost-1];
                    rho_w[j][size+ghost] = rho_w[j][size+ghost-1];
                    E[j][size+ghost] = E[j][size+ghost-1];
                }
                //Fix four conner 
                rho[0][0] = rho[1][0];
                rho[rho.length-1][0] = rho[rho.length-1][1];
                rho[0][rho.length-1] = rho[1][rho.length-1];
                rho[rho.length-1][rho.length-1] = rho[rho.length-2][rho.length-1];
                E[0][0] = E[1][0];
                E[rho.length-1][0] = E[rho.length-1][1];
                E[0][rho.length-1] = E[1][rho.length-1];
                E[rho.length-1][rho.length-1] = E[rho.length-2][rho.length-1];
                rho_u[0][0] = rho_u[1][0];
                rho_u[rho.length-1][0] = rho_u[rho.length-1][1];
                rho_u[0][rho.length-1] = rho_u[1][rho.length-1];
                rho_u[rho.length-1][rho.length-1] = rho_u[rho.length-2][rho.length-1];
                rho_v[0][0] = rho_v[1][0];
                rho_v[rho.length-1][0] = rho_v[rho.length-1][1];
                rho_v[0][rho.length-1] = rho_v[1][rho.length-1];
                rho_v[rho.length-1][rho.length-1] = rho_v[rho.length-2][rho.length-1];
                rho_w[0][0] = rho_w[1][0];
                rho_w[rho.length-1][0] = rho_w[rho.length-1][1];
                rho_w[0][rho.length-1] = rho_w[1][rho.length-1];
                rho_w[rho.length-1][rho.length-1] = rho_w[rho.length-2][rho.length-1];
                pressure[0][0] = pressure[1][0];
                pressure[rho.length-1][0] = pressure[rho.length-1][1];
                pressure[0][rho.length-1] = pressure[1][rho.length-1];
                pressure[rho.length-1][rho.length-1] = pressure[rho.length-2][rho.length-1];
                break;
            case 4: 
            	// All periodic
                for(int i = ghost;i <= size+ghost; i++) {
                    rho[0][i] = rho[size+ghost-1][i];
                    pressure[0][i] = pressure[size+ghost-1][i];
                    rho_u[0][i] = rho_u[size+ghost-1][i];
                    rho_v[0][i] = rho_v[size+ghost-1][i];
                    rho_w[0][i] = rho_w[size+ghost-1][i];
                    E[0][i] = E[size+ghost-1][i];
                    rho[size+ghost][i] = rho[ghost][i];
                    pressure[size+ghost][i] = pressure[ghost][i];
                    rho_u[size+ghost][i] = rho_u[ghost][i];
                    rho_v[size+ghost][i] = rho_v[ghost][i];
                    rho_w[size+ghost][i] = rho_w[ghost][i];
                    E[size+ghost][i] = E[ghost][i];
                }
                for(int i = ghost;i <= size+ghost; i++){
                    rho[i][0] = rho[i][size+ghost-1];
                    pressure[i][0] = pressure[i][size+ghost-1];
                    rho_u[i][0] = rho_u[i][size+ghost-1];
                    rho_v[i][0] = rho_v[i][size+ghost-1];
                    rho_w[i][0] = rho_w[i][size+ghost-1];
                    E[i][0] = E[i][size+ghost-1];
                    rho[i][size+ghost] = rho[i][ghost];
                    pressure[i][size+ghost] = pressure[i][ghost];
                    rho_u[i][size+ghost] = rho_u[i][ghost];
                    rho_v[i][size+ghost] = rho_v[i][ghost];
                    rho_w[i][size+ghost] = rho_w[i][ghost];
                    E[i][size+ghost] = E[i][ghost];
                }
                //Fix four conner 
                rho[0][0] = rho[1][0];
                rho[rho.length-1][0] = rho[rho.length-1][1];
                rho[0][rho.length-1] = rho[1][rho.length-1];
                rho[rho.length-1][rho.length-1] = rho[rho.length-2][rho.length-1];
                E[0][0] = E[1][0];
                E[rho.length-1][0] = E[rho.length-1][1];
                E[0][rho.length-1] = E[1][rho.length-1];
                E[rho.length-1][rho.length-1] = E[rho.length-2][rho.length-1];
                rho_u[0][0] = rho_u[1][0];
                rho_u[rho.length-1][0] = rho_u[rho.length-1][1];
                rho_u[0][rho.length-1] = rho_u[1][rho.length-1];
                rho_u[rho.length-1][rho.length-1] = rho_u[rho.length-2][rho.length-1];
                rho_v[0][0] = rho_v[1][0];
                rho_v[rho.length-1][0] = rho_v[rho.length-1][1];
                rho_v[0][rho.length-1] = rho_v[1][rho.length-1];
                rho_v[rho.length-1][rho.length-1] = rho_v[rho.length-2][rho.length-1];
                rho_w[0][0] = rho_w[1][0];
                rho_w[rho.length-1][0] = rho_w[rho.length-1][1];
                rho_w[0][rho.length-1] = rho_w[1][rho.length-1];
                rho_w[rho.length-1][rho.length-1] = rho_w[rho.length-2][rho.length-1];
                pressure[0][0] = pressure[1][0];
                pressure[rho.length-1][0] = pressure[rho.length-1][1];
                pressure[0][rho.length-1] = pressure[1][rho.length-1];
                pressure[rho.length-1][rho.length-1] = pressure[rho.length-2][rho.length-1];
                break;
        }
    }

    // HLLC flux calculation
    @Override
    public void updateFlux(){
        // Three signal speed
        double S_L, S_R, SS;
        // Centre of Box(Box frame), Edge of Box(Edge frame)
        // Box Frame: i and j locate the centre of target cell
        // Edge Frame: i and j represent i+1/2, j or i, j+1/2
        // inside a single cell, there are four flux
        for(int i = 0; i < size+ghost; i++){
            for(int j = 0; j < size+ghost; j++){
                // Here we calculate the flux at i, j+1/2 (y-flux) and i+1/2, j (x-flux)
                
                //Flux at i+1/2, j (x-flux)
                //Calculate the signal speed (Box frame to Edge frame)
                S_L = Math.min(
                        //Signal speed at \lambda-
                        rho_u[i][j]/rho[i][j]- Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_u[i+1][j]/rho[i+1][j]-Math.sqrt(gamma*pressure[i+1][j]/rho[i+1][j])
                );
                
                S_R = Math.max(
                        //Signal speed at \lambda-
                        rho_u[i][j]/rho[i][j] + Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_u[i+1][j]/rho[i+1][j] + Math.sqrt(gamma*pressure[i+1][j]/rho[i+1][j])
                );
                //Intermediate speed
                SS = ( pressure[i+1][j]-pressure[i][j] +
                    rho_u[i][j]*(S_L-rho_u[i][j]/rho[i][j]) -
                    rho_u[i+1][j]*(S_R-rho_u[i+1][j]/rho[i+1][j]) )/
                    ((S_L*rho[i][j]-rho_u[i][j])-(S_R*rho[i+1][j]-rho_u[i+1][j]));
                
                //HLLC flux determination
                //Here (i,j) [box frame] -> (i+1/2,j) [edge frame]
                if(0 <= S_L){
                    //Flux from i to i+1 through i+1/2
                    x_rhoF[i][j] = x_rho_Flux(i, j);
                    x_rho_uF[i][j] = x_rho_u_Flux(i, j);
                    x_rho_vF[i][j] = x_rho_v_Flux(i, j);
                    x_rho_wF[i][j] = x_rho_w_Flux(i, j);
                    x_EF[i][j] = x_E_Flux(i, j);
                }
                else if(0 >= S_R)
                {
                    //Flux from i+1 to i through i+1/2
                    x_rhoF[i][j] = x_rho_Flux(i+1, j);
                    x_rho_uF[i][j] = x_rho_u_Flux(i+1, j);
                    x_rho_vF[i][j] = x_rho_v_Flux(i+1, j);
                    x_rho_wF[i][j] = x_rho_w_Flux(i+1, j);
                    x_EF[i][j] = x_E_Flux(i+1, j);
                }else if(S_L <= 0 && 0 <= SS)
                {
                    //Flux from i intermediate to i+1 throught i+1/2
                    x_rhoF[i][j] = ( SS * ( S_L * rho[i][j] - x_rho_Flux(i, j) ) ) / ( S_L - SS );
                    x_rho_uF[i][j] = ( SS * ( S_L * rho_u[i][j] 
                            - x_rho_u_Flux(i, j) ) + S_L * 
                            (pressure[i][j] + rho[i][j] * 
                            (S_L-rho_u[i][j]/rho[i][j]) * (SS-rho_u[i][j]/rho[i][j]) )) 
                            / ( S_L - SS );
                    x_rho_vF[i][j] = ( SS * ( S_L * rho_v[i][j] - x_rho_v_Flux(i, j) )) 
                            / ( S_L - SS );
                    x_rho_wF[i][j] = ( SS * ( S_L * rho_w[i][j] - x_rho_w_Flux(i, j)) ) 
                            / ( S_L - SS );
                    x_EF[i][j] = ( SS * ( S_L * E[i][j] 
                            - x_E_Flux(i, j) ) + SS * S_L * 
                            (pressure[i][j] + rho[i][j] * 
                            (S_L-rho_u[i][j]/rho[i][j]) * (SS-rho_u[i][j]/rho[i][j]) )) 
                            / ( S_L - SS );
                }else if(SS <= 0 && 0 <= S_L)
                {
                    //Flux from i+1 intermediate to i throught i+1/2
                    x_rhoF[i][j] = ( SS * ( S_R * rho[i+1][j] - x_rho_Flux(i+1, j) ) ) / ( S_R - SS );
                    x_rho_uF[i][j] = ( SS * ( S_R * rho_u[i+1][j] 
                            - x_rho_u_Flux(i+1, j) ) + S_R * 
                            (pressure[i+1][j] + rho[i][j] * 
                            (S_R-rho_u[i+1][j]/rho[i+1][j]) * (SS-rho_u[i+1][j]/rho[i+1][j]) )) 
                            / ( S_R - SS );
                    x_rho_vF[i][j] = ( SS * ( S_R * rho_v[i+1][j] - x_rho_v_Flux(i+1, j)) ) 
                            / ( S_R - SS );
                    x_rho_wF[i][j] = ( SS * ( S_R * rho_w[i+1][j] - x_rho_w_Flux(i+1, j)) ) 
                            / ( S_R - SS );
                    x_EF[i][j] = ( SS * ( S_R * E[i+1][j] 
                            - x_E_Flux(i+1, j) ) + SS * S_R * 
                            (pressure[i+1][j] + rho[i][j] * 
                            (S_R-rho_u[i+1][j]/rho[i+1][j]) * (SS-rho_u[i+1][j]/rho[i+1][j]) )) 
                            / ( S_R - SS );
                }
                
                
                //Flux at i, j+1/2 (y-direction)
                //Calculate the signal speed (Box frame to Edge frame)
                S_L = Math.min(
                        //Signal speed at \lambda-
                        rho_v[i][j]/rho[i][j]- Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_v[i][j+1]/rho[i][j+1]-Math.sqrt(gamma*pressure[i][j+1]/rho[i][j+1])
                );
                
                S_R = Math.max(
                        //Signal speed at \lambda-
                        rho_v[i][j]/rho[i][j] + Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_v[i][j+1]/rho[i][j+1] + Math.sqrt(gamma*pressure[i][j+1]/rho[i][j+1])
                );
                //Intermediate speed
                SS = ( pressure[i][j+1]-pressure[i][j] +
                    rho_v[i][j]*(S_L-rho_v[i][j]/rho[i][j]) -
                    rho_v[i][j+1]*(S_R-rho_v[i][j+1]/rho[i][j+1]) )/
                    ((S_L*rho[i][j]-rho_v[i][j])-(S_R*rho[i][j+1]-rho_v[i+1][j+1]));
                
                //HLLC flux determination
                //Here (i,j) [box frame] -> (i+1/2,j) [edge frame]
                if(0 <= S_L){
                    //Flux from i to i+1 through i+1/2
                    y_rhoF[i][j] = y_rho_Flux(i, j);
                    y_rho_uF[i][j] = y_rho_u_Flux(i, j);
                    y_rho_vF[i][j] = y_rho_v_Flux(i, j);
                    y_rho_wF[i][j] = y_rho_w_Flux(i, j);
                    y_EF[i][j] = y_E_Flux(i, j);
                }
                else if(0 >= S_R)
                {
                    //Flux from i+1 to i through i+1/2
                    y_rhoF[i][j] = y_rho_Flux(i, j+1);
                    y_rho_uF[i][j] = y_rho_u_Flux(i, j+1);
                    y_rho_vF[i][j] = y_rho_v_Flux(i, j+1);
                    y_rho_wF[i][j] = y_rho_v_Flux(i, j+1);
                    y_EF[i][j] = y_E_Flux(i, j+1);
                }else if(S_L <= 0 && 0 <= SS)
                {
                    //Flux from i intermediate to i+1 throught i+1/2
                    y_rhoF[i][j] = ( SS * ( S_L * rho[i][j] - y_rho_Flux(i, j) ) ) / ( S_L - SS );
                    y_rho_uF[i][j] = ( SS * ( S_L * rho_u[i][j] - y_rho_u_Flux(i, j)) ) 
                            / ( S_L - SS );
                    y_rho_vF[i][j] = ( SS * ( S_L * rho_v[i][j] 
                            - y_rho_v_Flux(i, j) ) + S_L * 
                            (pressure[i][j] + rho[i][j] * 
                            (S_L-rho_v[i][j]/rho[i][j]) * (SS-rho_v[i][j]/rho[i][j]) )) 
                            / ( S_L - SS );
                    y_rho_wF[i][j] = ( SS * ( S_L * rho_w[i][j] - y_rho_v_Flux(i, j)) ) 
                            / ( S_L - SS );
                    y_EF[i][j] = ( SS * ( S_L * E[i][j] 
                            - y_E_Flux(i, j) ) + S_L * SS *
                            (pressure[i][j] + rho[i][j] * 
                            (S_L-rho_u[i][j]/rho[i][j]) * (SS-rho_u[i][j]/rho[i][j]) )) 
                            / ( S_L - SS );
                }else if(SS <= 0 && 0 <= S_L)
                {
                    //Flux from i+1 intermediate to i throught i+1/2
                    y_rhoF[i][j] = ( SS * ( S_R * rho[i][j+1] - y_rho_Flux(i, j+1) ) ) / ( S_R - SS );
                    y_rho_uF[i][j] = ( SS * ( S_R * rho_u[i][j+1] - y_rho_u_Flux(i, j+1)) ) 
                            / ( S_R - SS );
                    y_rho_vF[i][j] = ( SS * ( S_R * rho_v[i][j+1] 
                            - y_rho_v_Flux(i, j+1) ) + S_R * 
                            (pressure[i][j+1] + rho[i][j] * 
                            (S_R-rho_v[i][j+1]/rho[i][j+1]) * (SS-rho_v[i][j+1]/rho[i][j+1]) )) 
                            / ( S_R - SS );
                    y_rho_wF[i][j] = ( SS * ( S_R * rho_w[i][j+1] - y_rho_w_Flux(i, j+1)) ) 
                            / ( S_R - SS );
                    y_EF[i][j] = ( SS * ( S_R * E[i][j+1] 
                            - y_E_Flux(i, j+1) ) + S_R * SS *
                            (pressure[i][j+1] + rho[i][j] * 
                            (S_R-rho_v[i][j+1]/rho[i][j+1]) * (SS-rho_v[i][j+1]/rho[i][j+1]) )) 
                            / ( S_R - SS );
                }
            }
        }
    }

    public void updateHLLFlux(){
        // Three signal speed
        double S_L, S_R;
        // Centre of Box(Box frame), Edge of Box(Edge frame)
        // Box Frame: i and j locate the centre of target cell
        // Edge Frame: i and j represent i+1/2, j or i, j+1/2
        // inside a single cell, there are four flux
        for(int i = 0; i < size+ghost; i++){
            for(int j = 0; j < size+ghost; j++){
                // Here we calculate the flux at i, j+1/2 (y-flux) and i+1/2, j (x-flux)
                
                //Flux at i+1/2, j (x-flux)
                //Calculate the signal speed (Box frame to Edge frame)
                S_L = Math.min(
                        //Signal speed at \lambda-
                        rho_u[i][j]/rho[i][j]- Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_u[i+1][j]/rho[i+1][j]-Math.sqrt(gamma*pressure[i+1][j]/rho[i+1][j])
                );
                
                S_R = Math.max(
                        //Signal speed at \lambda-
                        rho_u[i][j]/rho[i][j] + Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_u[i+1][j]/rho[i+1][j] + Math.sqrt(gamma*pressure[i+1][j]/rho[i+1][j])
                );
                
                
                //HLL flux determination
                //Here (i,j) [box frame] -> (i+1/2,j) [edge frame]
                if(0 <= S_L){
                    //Flux from i to i+1 through i+1/2
                    x_rhoF[i][j] = x_rho_Flux(i, j);
                    x_rho_uF[i][j] = x_rho_u_Flux(i, j);
                    x_rho_vF[i][j] = x_rho_v_Flux(i, j);
                    x_rho_wF[i][j] = x_rho_w_Flux(i, j);
                    x_EF[i][j] = x_E_Flux(i, j);
                }
                else if(0 >= S_R)
                {
                    //Flux from i+1 to i through i+1/2
                    x_rhoF[i][j] = x_rho_Flux(i+1, j);
                    x_rho_uF[i][j] = x_rho_u_Flux(i+1, j);
                    x_rho_vF[i][j] = x_rho_v_Flux(i+1, j);
                    x_rho_wF[i][j] = x_rho_w_Flux(i+1, j);
                    x_EF[i][j] = x_E_Flux(i+1, j);
                }else 
                {
                    //Flux from i intermediate to i+1 throught i+1/2
                    x_rhoF[i][j] = ( S_R * x_rho_Flux(i, j) - 
                            S_L * x_rho_Flux(i+1, j)+
                            S_L*S_R*( rho[i+1][j] - rho[i][j] ) ) / ( S_R - S_L );
                    
                    x_rho_uF[i][j] = ( S_R * x_rho_u_Flux(i, j) - 
                            S_L * x_rho_u_Flux(i+1, j)+
                            S_L*S_R*( rho_u[i+1][j] - rho_u[i][j] ) ) / ( S_R - S_L );
                    x_rho_vF[i][j] = ( S_R * x_rho_v_Flux(i, j) - 
                            S_L * x_rho_v_Flux(i+1, j)+
                            S_L*S_R*( rho_v[i+1][j] - rho_v[i][j] ) ) / ( S_R - S_L );
                    x_rho_wF[i][j] = ( S_R * x_rho_w_Flux(i, j) - 
                            S_L * x_rho_w_Flux(i+1, j)+
                            S_L*S_R*( rho_w[i+1][j] - rho_w[i][j] ) ) / ( S_R - S_L );
                    x_EF[i][j] = ( S_R * x_E_Flux(i, j) - 
                            S_L * x_E_Flux(i+1, j)+
                            S_L*S_R*( E[i+1][j] - E[i][j] ) ) / ( S_R - S_L );
                }
                
               
                
                //Flux at i, j+1/2 (y-direction)
                //Calculate the signal speed (Box frame to Edge frame)
                S_L = Math.min(
                        //Signal speed at \lambda-
                        rho_v[i][j]/rho[i][j]- Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_v[i][j+1]/rho[i][j+1]-Math.sqrt(gamma*pressure[i][j+1]/rho[i][j+1])
                );
                
                S_R = Math.max(
                        //Signal speed at \lambda-
                        rho_v[i][j]/rho[i][j] + Math.sqrt(gamma*pressure[i][j]/rho[i][j]) , 
                        //Signal speed at \lambda+
                        rho_v[i][j+1]/rho[i][j+1] + Math.sqrt(gamma*pressure[i][j+1]/rho[i][j+1])
                );
                
                
                //Here (i,j) [box frame] -> (i+1/2,j) [edge frame]
                if(0 <= S_L){
                    //Flux from i to i+1 through i+1/2
                    y_rhoF[i][j] = y_rho_Flux(i, j);
                    y_rho_uF[i][j] = y_rho_u_Flux(i, j);
                    y_rho_vF[i][j] = y_rho_v_Flux(i, j);
                    y_rho_wF[i][j] = y_rho_w_Flux(i, j);
                    y_EF[i][j] = y_E_Flux(i, j);
                }
                else if(0 >= S_R)
                {
                    //Flux from i+1 to i through i+1/2
                    y_rhoF[i][j] = y_rho_Flux(i, j+1);
                    y_rho_uF[i][j] = y_rho_u_Flux(i, j+1);
                    y_rho_vF[i][j] = y_rho_v_Flux(i, j+1);
                    y_rho_wF[i][j] = y_rho_w_Flux(i, j+1);
                    y_EF[i][j] = y_E_Flux(i, j+1);
                }else 
                {
                    //Flux from i intermediate to i+1 throught i+1/2
                    y_rhoF[i][j] = ( S_R * y_rho_Flux(i, j) - S_L * y_rho_Flux(i, j+1)+
                            S_L*S_R*( rho[i][j+1] - rho[i][j] ) ) / ( S_R - S_L );
                    
                    y_rho_uF[i][j] = ( S_R * y_rho_u_Flux(i, j) - 
                            S_L* y_rho_u_Flux(i, j+1) +
                            S_L*S_R*( rho_u[i][j+1] - rho_u[i][j] ) ) / ( S_R - S_L );
                    y_rho_vF[i][j] = ( S_R * y_rho_v_Flux(i, j)- 
                            S_L* y_rho_v_Flux(i, j+1)+
                            S_L*S_R*( rho_v[i][j+1] - rho_v[i][j] ) ) / ( S_R - S_L );
                    y_rho_wF[i][j] = ( S_R * y_rho_w_Flux(i, j) - 
                            S_L * y_rho_w_Flux(i, j+1)+
                            S_L*S_R*( rho_w[i][j+1] - rho_w[i][j] ) ) / ( S_R - S_L );
                    y_EF[i][j] = ( S_R * y_E_Flux(i, j) - 
                            S_L * y_E_Flux(i, j+1)+
                            S_L*S_R*( E[i][j+1] - E[i][j] ) ) / ( S_R - S_L );
                }
            }
        }
    }
    
    // The following methods are used to find the states and flux
    private double u_squared(int i, int j){
        return (rho_u[i][j]*rho_u[i][j]+rho_v[i][j]*rho_v[i][j]+rho_w[i][j]*rho_w[i][j]) /
                (rho[i][j]*rho[i][j]);
    }
    
    private double E_finding(int i, int j){
        return rho[i][j] * u_squared(i,j) / 2 + pressure[i][j] / ( gamma - 1 ) ;
    }
    
    private double x_E_Flux(int i, int j){
        return rho_u[i][j] * ( E[i][j] + pressure[i][j] ) / rho[i][j];
    }
    
    private double y_E_Flux(int i, int j){
        return rho_v[i][j] * ( E[i][j] + pressure[i][j] ) / rho[i][j];
    }
    
    private double E_to_pressure(int i, int j){
        return ( E[i][j] - rho[i][j] * u_squared(i,j) / 2 ) * (gamma-1) ;
    }
    
    private double x_rho_Flux(int i, int j){
        return rho_u[i][j];
    }
    
    private double y_rho_Flux(int i, int j){
        return rho_v[i][j];
    }
    
    private double x_rho_u_Flux(int i, int j){
        return rho_u[i][j] * rho_u[i][j] / rho[i][j] + pressure[i][j];
    }
    
    private double y_rho_u_Flux(int i, int j){
        return rho_u[i][j] * rho_v[i][j] / rho[i][j];
    }
    
    private double x_rho_v_Flux(int i, int j){
        return rho_v[i][j] * rho_u[i][j] / rho[i][j];
    }
    
    private double y_rho_v_Flux(int i, int j){
        return rho_v[i][j] * rho_v[i][j] / rho[i][j] + pressure[i][j];
    }
    
    private double x_rho_w_Flux(int i, int j){
        return rho_w[i][j] * rho_u[i][j] / rho[i][j];
    }
    
    private double y_rho_w_Flux(int i, int j){
        return rho_w[i][j] * rho_v[i][j] / rho[i][j];
    }
    
    
    @Override
    public void updateState(){
        double rho_new[][] = new double[size+ghost*2][size+ghost*2];
        double rho_u_new[][] = new double[size+ghost*2][size+ghost*2];
        double rho_v_new[][] = new double[size+ghost*2][size+ghost*2];
        double rho_w_new[][] = new double[size+ghost*2][size+ghost*2];
        double pressure_new[][] = new double[size+ghost*2][size+ghost*2];
        double E_new[][] = new double[size+ghost*2][size+ghost*2];
        // i is x coordinate
        
        for(int i = ghost; i < size+ghost; i++){
            // j is y coordinate
            for(int j = ghost; j < size+ghost; j++){
                //Update the state with flux
                rho_new[i][j] = rho[i][j] - 
                        dt* (x_rhoF[i-ghost+1][j-ghost+1]-x_rhoF[i-ghost][j-ghost+1])/dx -
                        dt* (y_rhoF[i-ghost+1][j-ghost+1]-y_rhoF[i-ghost+1][j-ghost])/dy;              
                rho_u_new[i][j] = rho_u[i][j] - 
                        dt* (x_rho_uF[i-ghost+1][j-ghost+1]-x_rho_uF[i-ghost][j-ghost+1])/dx -
                        dt* (y_rho_uF[i-ghost+1][j-ghost+1]-y_rho_uF[i-ghost+1][j-ghost])/dy;
                rho_v_new[i][j] = rho_v[i][j] - 
                        dt* (x_rho_vF[i-ghost+1][j-ghost+1]-x_rho_vF[i-ghost][j-ghost+1])/dx -
                        dt* (y_rho_vF[i-ghost+1][j-ghost+1]-y_rho_vF[i-ghost+1][j-ghost])/dy;
                rho_w_new[i][j] = rho_w[i][j] - 
                        dt* (x_rho_wF[i-ghost+1][j-ghost+1]-x_rho_wF[i-ghost][j-ghost+1])/dx -
                        dt* (y_rho_wF[i-ghost+1][j-ghost+1]-y_rho_wF[i-ghost+1][j-ghost])/dy;
                E_new[i][j] = E[i][j] - 
                        dt* (x_EF[i-ghost+1][j-ghost+1]-x_EF[i-ghost][j-ghost+1])/dx -
                        dt* (y_EF[i-ghost+1][j-ghost+1]-y_EF[i-ghost+1][j-ghost])/dy;
                
                //Adiabatic process
                //pressure_new[i][j] = pressure[i][j]*Math.pow(rho_new[i][j]/rho[i][j], gamma);
                pressure_new[i][j] = E_to_pressure(i,j);
            }
        }
        
        rho = rho_new;
        rho_u = rho_u_new;
        rho_v = rho_v_new;
        rho_w = rho_w_new;
        E = E_new;
        pressure = pressure_new;
        
        setBoundary(3);
    }
    
    private double[] maxSpeed(){
        double maxSpeed[] = {0,0};
        for(int i = ghost; i < ghost+size; i++){
            for(int j = ghost; j < ghost+size; j++){
                double c = Math.sqrt(gamma*pressure[i][j]/rho[i][j]);
                if(maxSpeed[0] < c + Math.abs(rho_u[i][j]/rho[i][j]))
                    maxSpeed[0] = c + Math.abs(rho_u[i][j]/rho[i][j]);
                if(maxSpeed[1] < c+Math.abs(rho_v[i][j]/rho[i][j]))
                    maxSpeed[1] = c + Math.abs(rho_v[i][j]/rho[i][j]);
            }
        }
        if(maxSpeed[0] == 0 && maxSpeed[1] == 0){
            System.out.println("Error in calculating max Speed");
            System.exit(0);
        }
        return maxSpeed;
    }
    
    @Override
    public void checkCFL(){
        double[] maxSpeed = maxSpeed();
        dt = CFL * dx / (maxSpeed[0]+maxSpeed[1]);
    }

    public void iterate(double time) throws IOException{
        BufferedWriter bw = Files.newBufferedWriter(new File("xy2.txt").toPath());
        setIC();
        double t = 0;
        while(true){
            //Printing method
            /*if(t > time){
                bw.write(Double.toString(t));
                bw.newLine();
                for(int j = 0; j < size; j++){
                    bw.write(","+dx*j);
                }
                bw.newLine();
                for(int j = 0; j < size; j++){
                    bw.write(Double.toString(j*dy));
                    for(int k = 0; k < size; k++){
                        bw.write(","+rho[j+ghost][k+ghost]);
                    }
                    bw.newLine();      
                    bw.flush();
                }
                break;
            }*/
            
            // Printing method
            /*
            if(t > 0.2){
                double sum = 0;
                bw.write(Double.toString(t));
                bw.newLine();
                for(int j = 0; j < size; j++){
                    bw.write(dx*j+","+rho[4][j+ghost]);
                bw.newLine();
                }
                System.out.println(sum);
                break;
            }*/
            
            bw.write(Double.toString(t));
            bw.newLine();
            for(int j = 0; j < size; j++){
                bw.write(dx*j+","+rho_v[4][j+ghost]/rho[4][j+ghost]);
                bw.newLine();
            }
            
            if(t > time){
                bw.close();
                break;
            }
            
            updateHLLFlux();
            updateState();
            t += dt;
            System.out.println("iterating t = "+t);
            checkCFL();
        }
    }
}