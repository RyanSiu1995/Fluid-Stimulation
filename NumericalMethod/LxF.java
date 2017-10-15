package NumericalMethod;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class LxF {
    private double[] currentState, nextState;
    private final int size, ghost;
    private final double dx, dt;
    private double max, min;
    
    private class Flux{
        // a = dx/dt
        public double F1(double q0, double q1){
            return (q0*q0/2+q1*q1/2)/2-dx*(q1-q0)/(2*dt);
        }
        public double F2(double q_1, double q0){
            return (q_1*q_1/2+q0*q0/2)/2-dx*(q0-q_1)/(2*dt);
        }
        
        
        // a_{1/2} = max|f'(q)|
        public double LF1(double q0, double q1){
            return (q0*q0/2+q1*q1/2)/2-Math.max(Math.abs(q0), Math.abs(q1))*(q1-q0)/(2);
        }
        public double LF2(double q_1, double q0){
            return (q_1*q_1/2+q0*q0/2)/2-Math.max(Math.abs(q0), Math.abs(q_1))*(q0-q_1)/(2);
        }
    }
    
    private class MathTool{
        public int findLargest(double array[]){
            double largest = array[0];
            int largestIndex = 0;

            for(int i = 0; i < array.length; i++){
                if(array[i] > largest){
                    largest = array[i]; 
                    largestIndex =i;
                }  
            }
            return largestIndex;
        }
        
        public double sum_square(double[] q){
            double s = 0;
            for(int i = ghost; i < size+1; i++){
                s += q[i]*q[i];
            }
            return s;
        }
        
        public int findLastChange(double[] q, int start){
            int i;
            for(i = q.length-2; i > 0; i--){
                if(Math.abs(q[i+1] - q[i]) > 0.001){
                    break;
                }
            }
            return i;
        }
        
        public int find1stChange(double[] q, int start){
            int i;
            for(i = start; i < q.length; i++){
                if(q[i+1] - q[i] > 0.001){
                    break;
                }
            }
            return i;
        }
        
        public double sum(double[] q){
            double s = 0;
            for(int i = ghost; i < size+1; i++){
                s += q[i];
            }
            return s;
        }
        
    }
    
    public LxF(int size, double dx, double dt){
        this.size = size;
        this.ghost = 1; //1 ghost cell for LxF is enough
        this.dx = dx;
        this.dt = dt;
        currentState = new double[size+ghost*2];
        nextState = new double[size+ghost*2];
    }

    public void initializeICBC(char type, double max, double min){
        //initialize IC
        switch(type){
            case 's':
                for(int i = 0; i < size+1; i++){
                    if(i < 200)
                        currentState[i] = min;
                    else if(i < 400){
                        currentState[i] = min+(max-min)*Math.sin(((double)i-200)*Math.PI/200);
                    }
                    else
                        currentState[i] = min;
                }
        
                for(int i = 0; i < ghost; i++){
                    currentState[i] = min;
                    currentState[size+ghost+i] = min;
                }
                this.max = max;
                this.min = min;
                checkSetting(currentState);
                break;
            
            case 'S':
                for(int i = 0; i < size+1; i++){
                    if(i < 200)
                        currentState[i] = min;
                    else if(i < 400){
                        currentState[i] = max;
                    }
                    else
                        currentState[i] = min;
                }
        
                for(int i = 0; i < ghost; i++){
                    currentState[i] = min;
                    currentState[size+ghost+i] = min;
                }
                this.max = max;
                this.min = min;
                checkSetting(currentState);
                break;
        }
    }
    
    private void checkSetting(double[] currentState){
        if(checkEntropy()==false){
            System.out.println("Entropy condtion is not satisfied");
            System.out.println("Possible f' max ="+max);
            System.out.println("Possible f' min ="+min);
            System.out.println("Possible s ="+(max+min)/2);
            System.exit(0);
        }
        if(checkCFL(currentState)==false){
            MathTool m = new MathTool();
            double max_value = currentState[m.findLargest(currentState)];
            System.out.println("CFL condtion is not satisfied");
            System.out.println("f'(u)_max ="+max_value);
            System.out.println("dx/dt ="+dx/dt);
            System.exit(0);
        };
    }
    
    private boolean checkCFL(double[] currentState){
        MathTool m = new MathTool();
        double max_value = currentState[m.findLargest(currentState)];
        return (max_value <= dx/dt);
    }
    
    private boolean checkEntropy(){
        double s = (max - min)/2;
        return (max > s && min < s);
    }
    
    private double[] update(double[] q){
        Flux f = new Flux();
        double[] flux = new double[size+1];
        double[] next = q;
        
        for(int i = 0; i < size+1; i++){
            flux[i] = f.LF2(q[i], q[i+1]);
        }
        
        int j;
        for( j = ghost; j < size+ghost; j++){
                next[j] = q[j] - dt*(flux[j-ghost+1]-flux[j-ghost])/(dx);
                if(j == ghost)
                    next[size+ghost] = next[j];
                if(j == size)
                    next[0] = next[j];
        }
            //System.out.println(j);
           // this.nextState[0] = this.nextState[size];
            //this.nextState[size+ghost] = this.nextState[ghost];
       
        return next;
    }
    
    public void iterate() throws IOException, InterruptedException{
        BufferedWriter bw = Files.newBufferedWriter(new File("xy.txt").toPath());
        MathTool m = new MathTool();
        for(int i = 0; i < 1000; i++){
            /*bw.write("Particles");
            bw.newLine();
            for(int j = 0; j < size; j++){
                bw.write(j*dx+","+currentState[j+ghost]);
                bw.newLine();
            }*/
            
            /*w.write(Double.toString(dt*i)+","+Double.toString(t.MaxTracer(currentState, nextState)));
            bw.newLine();
            dummy = m.findLastChange(currentState,dummy);*/
            bw.write(Double.toString(m.sum_square(currentState)));
            bw.newLine();
            bw.flush();
            System.arraycopy(update(currentState), 0, currentState, 0, currentState.length);
            System.out.println("Finish iterating i = "+i);
            System.out.println(m.sum(currentState));
            
        }
        bw.close();
    }
}
