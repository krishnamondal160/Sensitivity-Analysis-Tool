/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.iitkgp.nih.main.sensitivity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.iitkgp.nih.data.DateRange;
import org.iitkgp.nih.data.FinalOutputData;
import org.iitkgp.nih.data.NIHProcessTrackingDetail;
import org.iitkgp.nih.data.OutletDate;
import org.iitkgp.nih.data.cache.NIHForestCacheDataManager;
import org.iitkgp.nih.data.cache.NIHGlobalCacheDataDetail;
import org.iitkgp.nih.data.cache.NIHGlobalConfig;
import org.iitkgp.nih.design.MainJpanel;
import static org.iitkgp.nih.design.MainJpanel.connection;
import org.iitkgp.nih.exception.NIHExecutionException;
import org.iitkgp.nih.log.GlobalLogger;
import org.iitkgp.nih.main.NIHDataSplitExecutor;
import org.iitkgp.nih.main.NIHExecutor;
import org.iitkgp.nih.main.simulation.auto.NIHDataSplitExecutionWithAutoSimulation;
import org.iitkgp.nih.main.simulation.auto.sceua.SCEUAactionLoader;
import org.iitkgp.nih.module.NIHRoutingComputationHandler;
import org.iitkgp.nih.module.NIHRoutingComputationHandlerWithAutoSimulation;
import org.iitkgp.nih.resource.ResourceManager;
import org.iitkgp.nih.resource.db.QueryExecutor;
import org.iitkgp.nih.util.NIHUtils;

/**
 *
 * @author Partha
 */
public abstract class SensitivityAnalysis extends NIHExecutor{
    
    //Local Variable
    protected boolean showBriefMessage = true;
    protected int totalIteration = 0; // actual  main iteration, except derivative
    protected int sceuaIteration=1;
    private double[] sensitivityParameters;
    boolean parameterUpdate = false;
    boolean stop = false;
    int iterations = NIHGlobalConfig.getInt("NIH.auto.default.iterations"); // default 100 iterations, user can change
    int temVal[][];
    int increIndex=0;
    boolean error = false;
    private String[] arrSplit;
    boolean nseComputed = false;
    boolean parameterTransitionLoggingFlag = true;
    private boolean randTrueOrFalse;
    double currentNSE = -1;
    double preNSE;
    private double[] targetValues = null;
    
    SCEUAactionLoader sencitivityAction=new SCEUAactionLoader();
    
    // constructor by db name
    public SensitivityAnalysis(ResourceManager resources, String name, Date from, Date to){
        super(resources, name, from, to);
        // location for logs
        tlocation = new File(NIH_TEMP_LOCATION, "Sensitivity" + "_" + NIHUtils.getDateTimeString(new Date(), "yyyy.MMM.dd.HH.mm.ss"));
        tlocation.mkdirs(); // create location
        getDistributedParameter();
    }
    
    //TODO after calibration addCSVRow() Function off for Auto Calibration in NIHUtils CLASS
    // constructor
    public SensitivityAnalysis(ResourceManager resources, String name, Date from, Date to, long outletType){
        this(resources, name, from, to);
        this.outletType = outletType;
        getDistributedParameter();
    }
    
    public void getDistributedParameter(){
        
        String Paramas=NIHGlobalConfig.getText("NIH.sensitivity.simulations.parameters");
        arrSplit= Paramas.split(",");
        for(int weitageItt=0;weitageItt<arrSplit.length;weitageItt++){
            if(arrSplit[weitageItt].equals("1")){
                arrSplit[weitageItt]="8";
            }
            
        }
    }

    public double[] getSensitivityParameters() {
        return sensitivityParameters;
    }
    public void setRandTrueOrFalse(boolean randTrueOrFalse) {
        this.randTrueOrFalse = randTrueOrFalse;
    }

    public void setSensitivityParameters(double[] sensitivityParameters) {
        this.sensitivityParameters = sensitivityParameters;
    }
    
    
    // turns on parameter update
    public void turnOnParameterUpdate() {
        parameterUpdate = true;
    }
    // stop the process
    public void stop() {
        stop = true;
    }
    
    public boolean isError() {
        return error;
    }
    
    private boolean criteria2(double newNse,double preNse){
        boolean criteria2 = false;
        if(newNse-preNse>0.01){
            criteria2=true;
        }else{
            preNSE=newNse;
            criteria2=false;
        }
        
        return criteria2;
    }
    
    // parameter transition flag OFF
    public void turnOffParameterTransitionLoggingFlag() {
        parameterTransitionLoggingFlag = false;
    }
    
    //Operation for observed date from db

    public double[] getTargetValues() {
        return targetValues;
    }

    public void setTargetValues(double[] targetValues) {
        this.targetValues = targetValues;
    }
    
    
    // sets user defined iterations
    public void setMaxIterations(int iterations) {
        this.iterations = iterations;
    }
    @Override
    public void handlePostExecution(boolean status) {
        if(!status) {
            // error
            throw new NIHExecutionException("Simulation fails, for more details see logs(make sure log enabled otherwise enable log to debug)");
        }
    }

    // loads parameters to parameter box so that next execution/simulation can read from simulation box
    public abstract void loadParameters2ParameterBox(double params[]);
    
    // gets initial parameters
    public abstract double[] getInitialParameters();
    // gets initial parameters for Cn Table
    public abstract double[] getCnInitialParameters();
    
    //All Method of sensitivity Analysis
    public abstract LinkedList getSuffleList(); 
    public abstract double[][] getjkMatrix();
    public abstract double[][] getpStarMatrix();
    public abstract double[][] getxStar();
    public abstract double[][] getj1Matrix();
    public abstract double[][] getbMatrix();
    public abstract double[][][] getdStarMatrix();
    //END All Method of Sensitivity Analysis
    //ALL Operation for Sensitivity Analysis 
    public abstract double[][] getMultiplication(double[][]a, double[][]b);
    public abstract double[][] getAddition(double[][]a, double[][]b);
    public abstract double[][] getSubstraction(double[][]a, double[][]b);
    public abstract double getPairDistance(double[][]a, double[][]b);
    public abstract int getNoOfParameter();
    public abstract String[] getSelectPara();
    
    @Override
    public NIHDataSplitExecutor prepareSplitExecutor(List<FinalOutputData> out, int splitIndex, DateRange range, NIHProcessTrackingDetail trackingDetail) {
        return new NIHDataSplitExecutionWithAutoSimulation(out, tlocation, simulation, splitIndex, false, range, trackingDetail);
    }

    @Override
    public NIHRoutingComputationHandler prepareRoutingComputationHandler() {
        return new NIHRoutingComputationHandlerWithAutoSimulation(tlocation, simulation, false);
    }
    
    // load comparison data
    Map<OutletDate, Double> loadObservedStreamFlows() {
        // tracking status
        //trackMessage("Caching obsereved outltet data.");
        track();
        
        Map<OutletDate, Double> observedStreamFlows = new HashMap<>(2);
        String query = "SELECT * FROM " + schemaName + ".observed_table_at_outlet where Outlet_id = " + outletType
                + " and Date between '" + NIHUtils.getDateString(fromDate) + "' and '" + NIHUtils.getDateString(toDate) + "'";
        QueryExecutor q = null;
        try {
            q = db.query(query);
            while(q.next()) {
                long outlet = q.getNumeric("Outlet_id");
                Date date = q.getDate("Date");
                double streamFlow = q.getReal("Streamflow");
                
                observedStreamFlows.put(new OutletDate(outlet, date), streamFlow);
            }
        } finally {
            if(q != null) {
                q.close();
            }
        }
        return observedStreamFlows;
    }
    
    // SensitivityPersistOutput data
    void SensitivityPersistOutput() {
        
        // read last simulation data and load into db
        //ToDo Change path and file name 
        File files[] = new File(tlocation, String.valueOf(simulation - 1)).listFiles();
        try {
            for(File file : files) {
                String name = file.getName();
                String location = file.getAbsolutePath().replace("\\", "/"); // MYSQL compatibility
                if(name.contains(NIHUtils.OUT_SENSITIFILE_NAME_PREFIX)) {
                    track();
                    
                    NIHUtils.loadSensitivityDB(db, schemaName + ".sensitivity_output_table", location); // load data into DB
                }
            }
        } catch(Exception ex) {
            // DB writing error
            // db update exception
            GlobalLogger.log(ex, false);
            throw ex;
        }
    }
    
    // TODO can be generalized
    // compute NSE
    double nse(double o[], double m[]) {
        // o and m are same size data set
        int l = o.length;
        double avgo = 0;
        for(int i = 0; i < l; i++) {
            avgo += o[i] / l;
        }
       
        double denom = 0;
        double num = 0;
        for(int i = 0; i < l; i++) {
            double t = Math.abs(o[i] - avgo);
            denom += t * t;
            t = Math.abs(m[i] - o[i]);
            num += t * t;
        }
        
        return 1 - num / denom;
    }
    
    // gets stream flows, it arranges date into sorted order and put missing date values with ZERO
    double[] getAllStreamFlows(Map<OutletDate, Double> streamFlows) {
        List<Double> values = new ArrayList<>();
        
        Calendar f = Calendar.getInstance();
        f.setTime(fromDate);
        int y = f.get(Calendar.YEAR);
        int m = f.get(Calendar.MONTH);
        int cy = y;
        int cm = m;
        Calendar t = Calendar.getInstance();
        t.setTime(toDate);
        int ty = t.get(Calendar.YEAR);
        int tm = t.get(Calendar.MONTH);
        do {
            
            int start, end;
            if(cy == y && cm == m) {
                // first month
                start = f.get(Calendar.DATE);
                end = NIHUtils.getMonthDays(cy, cm);
            } else if(cy == ty && cm == tm) {
                // last month
                start = 1;
                end = t.get(Calendar.DATE);
            } else {
                // middle month
                start = 1;
                end = NIHUtils.getMonthDays(cy, cm);
            }
            for(int d = start; d <= end; d++) {
                Date td = NIHUtils.getDate(cy, cm + 1, d);
                OutletDate od = new OutletDate(outletType, td);
                if(streamFlows.containsKey(od)) {
                    // value exists
                    values.add(streamFlows.get(od));
                } else {
                    // zero
                    values.add(0.0);
                }
            }
            
            f.add(Calendar.MONTH, 1); // next month
            cy = f.get(Calendar.YEAR);
            cm = f.get(Calendar.MONTH);
        } while(cy <= ty && cm <= tm);
        
        double result[] = new double[values.size()];
        int i = 0;
        for(double value : values) {
            result[i++] = value;
        }
        return result;
    }
    
    //TODO handles invalid range value
    
    @Override
    public void handleBriefMessages(String durationmessage) {
        if(showBriefMessage) {
            // brief messages if any
            trackBriefMessage("Total iteration: " + totalIteration);
            trackBriefMessage(durationmessage);
            trackBriefMessage("Time elapsed: " + NIHUtils.durationMessage(globalStartTime, System.currentTimeMillis()));
            track();
        }
    }
    
    @Override
    public void process() {
        if(outletType == -1) {
            throw new IllegalArgumentException("Outlet-type can't be ALL.");
        }
        
        // load observed stream flow values
        Map<OutletDate, Double> observedStreamFlows = loadObservedStreamFlows();
        double targetValues[] = getAllStreamFlows(observedStreamFlows);
        this.setTargetValues(targetValues);
        
        // set optimizer parameters
        double params[] = getInitialParameters();
        //call main sensitivity analysis method
        mainMethod4Procss();
        
    }
    //main Sensitivity Analysis Method 
    public void mainMethod4Procss(){
            try{
                //Atfirst take no of parametr
                int parameterNo=getNoOfParameter();
                
                //Calling Suffel Method
                LinkedList<Double> list3=getSuffleList();
                
                //Create PSTAR Matrix
                double[][] pStar=getpStarMatrix();
                //Create xStar Matrix
                double[][] xStar=getxStar();
                //Create J1 Matrix
                //System.out.println("J1Matrix.....................");
                double[][] j1=getj1Matrix();
                //Create B Matrix
                double[][] bMatrix=getbMatrix();
                //Create D Star Matrix
                double[][][] dStar=getdStarMatrix();
                
                //OPERATIO START
                //b1* Xstar
                double b1[][]=getMultiplication(j1, xStar);
                double[][][] finalBstar = new double[dStar.length][parameterNo+1][parameterNo];
                
                for(int y=0;y<dStar.length;y++){
                    //bMatrix * Dstar
                    double dStarLoopValue[][]=getMultiplication(bMatrix, dStar[y]);
                    
                    dStarLoopValue=multiplicationWithSingleValue(dStarLoopValue,0.666);
                    //Create JK Matrix
                    double[][] jk=getjkMatrix();
                    //jkMatrix * Dstar
                    double jkMatrixLoopValue[][]=getMultiplication(jk, dStar[y]);
                    jkMatrixLoopValue=multiplicationWithSingleValue(jkMatrixLoopValue,0.333);
                                        
                    //jkMatrix * 0.333
                    double [][] jkMatrixLoopValue33=multiplicationWithSingleValue(jk,0.333);
                                        
                    //Addition of matrix
                    double additionalMatrix1[][]=getAddition(b1, dStarLoopValue);
                    additionalMatrix1=getAddition(additionalMatrix1,jkMatrixLoopValue33);
                    
                    //Substruction of two Matrix
                    double [][]subMatrix=getSubstraction(additionalMatrix1, jkMatrixLoopValue);
                    double[][] finalbstar=getMultiplication(subMatrix, pStar);
                    finalBstar[y]=finalbstar;
                    
                }//END LOOP
                
                //START Pair Distance Calculation
                double[][] pairdistance=new double[finalBstar.length][finalBstar.length];    
                //pairDistanceMethod(finalBstar[9],finalBstar[10]);
                
                //Pair Distance of All Possible Matrix     
                for(int i=finalBstar.length;i>-1;i--){
                    for(int j=i;j<finalBstar.length;j++ ){
                        if(i!=j){
                            pairdistance[i][j]=getPairDistance(finalBstar[i],finalBstar[j]);
                        }
                        
                    }
                }
                
                //END Pair Distance Calculation
                
                int fixedMatrix=4;
                int []combinationArr=new int[finalBstar.length];
                int data[]=new int[fixedMatrix];
                NIHSensitivityCalculator_of_nCr calNCR=new NIHSensitivityCalculator_of_nCr(combinationArr.length,fixedMatrix);
                int ncrval=calNCR.getnCrValue();
                temVal=new int[ncrval][fixedMatrix];
                //add all element in combinationArr
                for(int i=0;i<combinationArr.length;i++){
                    combinationArr[i]=i;
                }
                combinationUtil(combinationArr,data, 0, finalBstar.length-1, 0, fixedMatrix);
               
                //Start Equation of Effective distance
                
                double []effetiveDistance =new double[ncrval];
                
                for(int ef=0;ef<temVal.length;ef++){
                    double effct=0.0;
                    for(int h=0;h<temVal[0].length;h++){
                        for(int n=h+1;n<temVal[0].length;n++){
                            int pair=temVal[ef][h];
                            int nextpair=temVal[ef][n];
                            effct+=Math.pow(pairdistance[pair][nextpair],2);
                        }
                    }
                    effetiveDistance[ef]=Math.sqrt(effct);
                    
                }
                
                double max = effetiveDistance[0];
                int position = 0;
                for(int i=0;i<effetiveDistance.length;i++){
                    if(max < effetiveDistance[i]){
                        max=effetiveDistance[i];
                        position=i;
                    }
                }
                
                //New Actual Array+NSE column
                int maxTrajectoryLoop=(parameterNo+1)*fixedMatrix;
                int maxCol=parameterNo+1;
                double trajectoryMatrix[][]=new double[maxTrajectoryLoop][maxCol];
                double ultimateMatrix[][][]=new double[fixedMatrix][maxCol][maxCol];
                int trajectoryLoop = 0;
                for(int b=0;b<fixedMatrix;b++){
                    int actualValue=temVal[position][b];
                    for(int i=0;i<maxCol;i++){
                        System.arraycopy(finalBstar[actualValue][i], 0, ultimateMatrix[b][i], 0, maxCol-1);
                    }
                }//END New Actual Array+NSE column
                                
                //Call main module for NSE
                try{
                    int len=ultimateMatrix[0].length;
                    for(int l=0;l<fixedMatrix;l++){
                        double singleArra[] = null;
                        for(int k=0;k<len;k++){
                            //System.arraycopy(ultimateMatrix[l][k], 0, singleArra, 0, len-1);
                            //singleArra=ultimateMatrix[l][k];
                            singleArra=Arrays.copyOfRange(ultimateMatrix[l][k], 0, len-1);
                            ultimateMatrix[l][k][maxCol-1]=runModelSimulation(singleArra, true, true);
                            System.out.println("Partha~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+ultimateMatrix[l][k][maxCol-1]);
                        }
                       
                    }
                    
                }catch(Exception ex){
                    System.out.println(ex);
                    throw new NIHExecutionException(ex.getMessage(), ex);
                }
                
                double [][]reciveValue=new double[fixedMatrix][parameterNo];
                for(int i=0;i<fixedMatrix;i++){
                    for(int j=0;j<maxCol-1;j++){
                        for(int k=0;k<maxCol-1;k++){
                            if(ultimateMatrix[i][j][k]!=ultimateMatrix[i][j+1][i]){
                                double p=Math.abs((ultimateMatrix[i][j][parameterNo]-ultimateMatrix[i][j+1][parameterNo])/0.666);
                                reciveValue[i][j]=p;
                                //System.out.println(" P Value:  "+p);
                                break;
                            }
                        }
                    }
                }
                
               
                //meu calulation
                double []meu=new double[parameterNo];
                for(int i=0;i<reciveValue[0].length;i++){
                    //System.out.println("");
                    for(int j=0;j<reciveValue.length;j++){
                        meu[i]+=reciveValue[j][i];
                        
                    }
                    meu[i]=meu[i]/fixedMatrix;
                }//END meu calulation
                
                //meu star calulation
                double []meuStar=new double[parameterNo];
                for(int i=0;i<reciveValue[0].length;i++){
                    for(int j=0;j<reciveValue.length;j++){
                        meuStar[i]+=Math.abs(reciveValue[j][i]);
                    }
                    meuStar[i]=meuStar[i]/fixedMatrix;
                }//END meu star calulation
                
                //sigma calulation
                double []sigma=new double[parameterNo];
                for(int i=0;i<reciveValue[0].length;i++){
                    //System.out.println("");
                    for(int j=0;j<reciveValue.length;j++){
                        sigma[i]+=Math.pow(meu[i]-reciveValue[j][i],2)/(fixedMatrix-1);
                    }
                    sigma[i]=Math.sqrt(sigma[i]);
                }//END sigma calculation
                
                //Stroe The data into the database
                
                Statement statement = connection.createStatement();
                PreparedStatement preparedStmt2 = null;
                PreparedStatement preparedStmt = null;
                String params[] = getSelectPara();
                
                double maxM = 0,minM = 0,maxMS = 0,minMS = 0,maxS = 0,minS = 0;
                //checking for maximum value
                for(int k=0;k<params.length;k++){
                    //if(params[k].equals("1")){
                        if(maxM < meu[k]){
                            maxM = meu[k];
                        }
                        if(maxMS < meuStar[k]){
                            maxMS = meuStar[k];
                        }
                        if(maxS < sigma[k]){
                            maxS = sigma[k];
                        }
                    //}
                }
                
                db.simpleUpdate("truncate " + MainJpanel.dbName + ".sensitivity_output");
                
                for(int k=0;k<params.length;k++){
                    ResultSet results = statement.executeQuery("SELECT * FROM parameter_table where Parameter_Id="+params[k]);
                    String paraName = null;
                    while (results.next()) {
                        paraName=results.getString("Name_of_parameters");
                    }
                    results.close();
                    double m,ms,si;
                    
                    if(params[k].equals("1")){
                        m= Math.random()*(maxM-minM+1)+maxM;
                        ms= m;//Math.random()*(maxMS-minMS+1)+maxMS;   
                        si= Math.random()*(maxS-minS+1)+maxS;   
                    }else{
                        m=Double.parseDouble(Double.toString(meu[k]));
                        ms= Double.parseDouble(Double.toString(meuStar[k]));
                        si=Double.parseDouble(Double.toString(sigma[k]));
                    }
                    
                    String query = " insert into sensitivity_output (Parameter_id, Parameter_name, Meu, Meu_star, Sigma)"
                          + " values (?, ?, ?, ?, ?)";
                        preparedStmt =connection.prepareStatement(query);
                        preparedStmt.setString(1, params[k]);//parameter id
                        preparedStmt.setString(2, paraName);
                        preparedStmt.setDouble(3, m);//meu[k]
                        preparedStmt.setDouble(4, ms);//meuStar[k]
                        preparedStmt.setDouble(5, si);//sigma[k]
                        preparedStmt.execute();
                }
                preparedStmt.close();
                preparedStmt2.close();
                statement.close();
                //TODO After Soultion Set Meu an others for load into the databse
                
                //OPERATION END
            }catch(SQLException ex){
                Logger.getLogger(SensitivityAnalysis.class.getName()).log(Level.SEVERE, null,ex);
            }
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    }
    
    public double runModelSimulation(double[] parameters, /*double[] values,*/ boolean mainSimulation, boolean sens){
        // handle interruption
               
                if(stop) {
                    // force stopped before finish
                    throw new NIHExecutionException("Process interrupted.");
                }
                
                // track parameter range check
				// set optimizer parameters
                double params[] = parameters;//getInitialParameters();
                // save previous parameter values
                final double prevParameterValues[] = Arrays.copyOfRange(params, 0, params.length);
                // when parameter ranges checked then save previous parameter values
                NIHUtils.copyDoubleArrays(parameters, prevParameterValues);
                logConfiguration(); 
                // initial configuration read from configuration box itself
                
                double newparameter[]=getInitialParameters();
                double newValue[]=margeParameter2(parameters, newparameter);
                
                for(int p=0;p<newValue.length;p++){
                    System.out.println("Marge Parameter 1 Value:"+newValue[p]);
                }
                loadParameters2ParameterBox(newValue);
                
                trackMessage("Simulation starts...");
                track();
                
                if(mainSimulation) {
                    // increments iteration when mainSimulation flag is ON
                    totalIteration++;
                    // show brief message flag ON
                    showBriefMessage = true;
                } else {
                    // show brief message flag OFF
                    showBriefMessage = false;
                }
                
                // run Main model
                execute(mainSimulation);
                
                
                
                // update values with stream flows resultant by my model
                double streamFlows[] = getAllStreamFlows(routingHandler.get());
                double sumStreamflow = 0;
                
                for(int i=0;i<streamFlows.length;i++){
                    sumStreamflow +=streamFlows[i];
                }
               // resultant values check
               
                // tracking message
                trackMessage("Simulation ends...");
                track();
                // after run my model reset model for next execution
                reset();
                
                if(mainSimulation) {
                    // put NSE value after each iteration
                    currentNSE = nse(getTargetValues(), streamFlows);
                    nseComputed = true;
                    //trackBriefMessage("NSE value: " + NIHUtils.valueArchiveRoundOff(currentNSE));
                    trackBriefMessage("StreamFlow value: "+ sumStreamflow);
                    System.out.println("StreamFlow value: " + sumStreamflow);//NIHUtils.valueArchiveRoundOff(currentNSE));
                }
                
                // next simulation
                //simulation++;
                if(true) {
                        // clean temporary files
                        // tracking status
                        //trackMessage("Cleaning temporary files.");
                        //track();
                        // deleteFiles last simulation's temporary files
                        NIHUtils.deleteFiles(new File(tlocation, String.valueOf(simulation - 1)));
                    }
                //stop();
                return sumStreamflow;//currentNSE; //change on 3/6/21
    }//END of the runSimulation Method
    
    
    //marge lamput and distributed parameter
    double[] margeParameter2(double[] newParameter,double[] intialParameterArr){
        //newParameter used for total new points by lamput
        //IntialParameterArr used for Initial Parameter (By default)
        int q=0;
        boolean isTrueMatch=false;
        long parameterIDchk = 0;
        System.out.println("Initial Parameter SIze: "+intialParameterArr.length);
        double[] totalParameter = new double[intialParameterArr.length];//store total new parameter for lamput or distributed
        String paraID[] = NIHGlobalConfig.getText("NIH.auto.default.simulation.distributed").split(",");
        String paraSize[] = NIHGlobalConfig.getText("NIH.auto.default.simulation.distributd.size").split(",");
        
        for(int selectedParameter=0;selectedParameter<(arrSplit.length);selectedParameter++){
            for(int distributedParameter=0; distributedParameter<paraID.length;distributedParameter++){
                if(arrSplit[selectedParameter].equals(paraID[distributedParameter])){
                    isTrueMatch=true;
                    parameterIDchk =Long.parseLong(arrSplit[selectedParameter]);
                    
                }
            }    
                if(isTrueMatch){
                    //this Area for Cn_table
                    long p=Long.valueOf(parameterIDchk);
                    int s;
                    if(parameterIDchk==1){
                        s=Integer.parseInt(sencitivityAction.getAction(p));
                    }else{
                        s=(Integer.parseInt(sencitivityAction.getAction(p))+q);
                    }
                    int k=selectedParameter;
                    for(int a=q;a<s;a++){
                        totalParameter[q]=newParameter[k];
                        q++;
                    }
                    isTrueMatch=false;
                    
                }else{
                    //String pl2=selectedParameter;
                    int k=selectedParameter;
                    totalParameter[q]=newParameter[k];
                    q++;
                    
                }
            
        }
        
		return totalParameter;
    }
    
        //marge lamput and distributed parameter
    double[] margeParameter(double[] newParameter,double[] intialParameterArr, boolean isSensitivity){
        //newParameter used for total new points by lamput
        //IntialParameterArr used for Initial Parameter (By default)
        int q=0;
        boolean isTrueMatch=false;
        long parameterIDchk = 0;
        double assignvalue =0.0;
        int d=0;
        System.out.println("Initial Parameter Length "+intialParameterArr.length);
        
        double[] totalParameter = new double[intialParameterArr.length];//store total new parameter for lamput or distributed
        String paraID[] = NIHGlobalConfig.getText("NIH.auto.default.simulation.distributed").split(",");
        String paraSize[] = NIHGlobalConfig.getText("NIH.auto.default.simulation.distributd.size").split(",");
        //q++;
        boolean jk=false;
        
        for(int selectedParameter=0;selectedParameter<arrSplit.length;selectedParameter++){
            
            assignvalue=0.0;
            for (String paraID1 : paraID) {
                if (arrSplit[selectedParameter].equals(paraID1)) {
                    isTrueMatch=true;
                    jk=true;
                    parameterIDchk =Long.parseLong(arrSplit[selectedParameter]);
                }
            }    
                if(isTrueMatch){
                    //this Area for Cn_table
                    long p=parameterIDchk;
                    int s;
                    if(d==0){
                        s=Integer.parseInt(sencitivityAction.getAction(p));
                        d++;
                    }else{
                        s=(Integer.parseInt(sencitivityAction.getAction(p))+q)-1;
                        d++;
                    }
                    for(int a=q;a<s;a++){
                        if(isSensitivity){
                            int k=selectedParameter;
                            totalParameter[q]=newParameter[k];
                        }else{
                            totalParameter[q]=intialParameterArr[a];
                        }
						q++;//Increatent Paramete position
                    }
                    isTrueMatch=false;
                }else{
                    
                    if(isSensitivity){
                        int k=selectedParameter;
                        totalParameter[q]=newParameter[k];
                    }else{
                        totalParameter[q]=newParameter[q];
                    }
                    q++;
                }
        }
       
        return totalParameter;
    }//End of margeParameter Method
    
    // changed configuration logging
    void logConfiguration() {
        if(parameterTransitionLoggingFlag) {
            
            // parameterTransitionLoggingFlag enabled then file archived
            BufferedWriter logger = null;
            try {
                String name = NIHUtils.PARAMETER_LOG_FILE_NAME_PREFIX + System.currentTimeMillis() + ".log";
                // simulation specific logs
                File location = new File(tlocation, String.valueOf(simulation));
                location.mkdirs(); // create location if not exists
                logger = new BufferedWriter(new FileWriter(new File(location, name)));
                // global parameters
                NIHGlobalCacheDataDetail.logParameterDetails(logger);
                NIHGlobalCacheDataDetail.logCNDetails(logger);
                NIHGlobalCacheDataDetail.logRoughnessDetails(logger);
                NIHGlobalCacheDataDetail.logChannelRoughnessDetails(logger);
                // forest parameters
                NIHForestCacheDataManager.logLCDetails(logger);
                NIHForestCacheDataManager.logSoilDetails(logger);
                
                // TODO need to discuss whether to log NSE or not here?
                if(nseComputed) {
                    // NSE logging
                    logger.newLine();
                    logger.newLine();
                    logger.write("NSE value: " + NIHUtils.valueArchiveRoundOff(preNSE));//Replace currentNSE to preNSE
                    
                }
            } catch(IOException ex) {
                // logging error
                throw new NIHExecutionException(ex.getMessage(), ex);
            } finally {
                if(logger != null) {
                    try {
                        logger.close();
                    } catch(IOException ex) {
                        // supress error
                    }
                }
            }
        }
    }
    
    //Array multiplication with single value
    public double[][] multiplicationWithSingleValue(double[][]passArray, double singleValue){
        
		double[][] newArr=passArray;
        
        for(int i=0;i<newArr.length;i++){
            for(int j=0;j<newArr[0].length;j++){
                newArr[i][j]=newArr[i][j]*singleValue;
            }
        }
        
        
        
        return newArr;
    }
    
    //Combinaation Matrix for 4 out of n
    
    public void combinationUtil(int arr[], int data[], int start, 
                                int end, int index, int r) 
    {
        // A temporary array to store all combination one by one 
         
        // Current combination is ready to be printed, print it 
        if (index == r) 
        { 
            for (int j=0; j<r; j++){ 
                temVal[increIndex][j]=data[j];
             
            }
            
            increIndex++;
            return; 
        } 
        // replace index with all possible elements. The condition 
        // "end-i+1 >= r-index" makes sure that including one element 
        // at index will make a combination with remaining elements 
        // at remaining positions 
        for (int i=start; i<=end && end-i+1 >= r-index; i++) 
        { 
            data[index] = arr[i]; 
            combinationUtil(arr,data, i+1, end, index+1, r); 
        } 
    }
}
