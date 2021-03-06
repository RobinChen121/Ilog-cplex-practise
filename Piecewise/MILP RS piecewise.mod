/*********************************************
 * OPL 12.6.0.0 Model
 * Author: ltopuser
 * Creation Date: 16 Oct 2015 at 10:44:43
 *********************************************/
//parameters
int nbmonths=...;
range months=1..nbmonths;
float fc=...;
float h=...;
float p=...;
float meandemand[months]=...;
float std_demand[months]=...;
float initialstock=...;

int nbpartitions=...;
range partitions=1..nbpartitions;
float means[partitions]=...;
float prob[partitions]=...;
float error=...;

//variables
dvar float stock[0..nbmonths];
dvar float+ stockhlb[0..nbmonths];
dvar float+ stockplb[0..nbmonths];
dvar boolean purchase[months];
dvar boolean P[months][months];

//float mean_matrix[i in months, j in months] = sum(m in i..j) meandemand[m];
float std_matrix[i in months, j in months] = sqrt(sum(m in i..j) pow(std_demand[m],2));

//objective function
minimize sum(t in months)( fc*purchase[t]+h*stockhlb[t]+p*stockplb[t]);

//constraints
subject to{
 stock[0]==initialstock;
 stockhlb[0]==maxl(stock[0],0);
 stockplb[0]==maxl(-stock[0],0);
 
 //purchase[1]==0; // x[1] can be zero

 forall(t in months) 
   		purchase[t] == 0 => stock[t]+meandemand[t]-stock[t-1] == 0;     
   		   
forall(t in months)
  		stock[t]+meandemand[t]-stock[t-1]>=0;
  
forall (t in months)
  		sum(j in 1..t)P[j][t]==1;
 
forall (t in months, j in 1..t)
  		P[j][t]>=purchase[j]-sum(k in j+1..t)purchase[k];
  
forall (t in months)
   		sum(k in 1..t) purchase[k] == 0 => P[1][t] == 1;  
  
/**
* Original formulation as in (Rossi et al., 2015)    
*/ 
forall(t in months, p in partitions)
  stockhlb[t]>=sum(k in 1..p)prob[k]*stock[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P[j][t]);

  forall(t in months) stockhlb[t] >= (sum(j in 1..t) error*std_matrix[j][t]*P[j][t]);

/**
* Piecewise-based formulation (requires an even number of partitions: ftoi(round(nbpartitions/2)) )
* 
* inside is standard normal distribution
* S = stock[t] + meand, S - meand= stock[t]
*/
forall(t in months, j in 1..t) 
  P[j][t] == 1 => stockhlb[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {((sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 1} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock[t]/std_matrix[j][t]);
    
/**
* Original formulation as in (Rossi et al., 2015)    
*/ 
forall(t in months, p in partitions)
  stockplb[t]>=-stock[t]+sum(k in 1..p)prob[k]*stock[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P[j][t]);
 
  forall(t in months) stockplb[t] >= - stock[t] + (sum(j in 1..t) error*std_matrix[j][t]*P[j][t]);

 
/**
* Piecewise-based formulation (requires an even number of partitions: ftoi(round(nbpartitions/2)) )
* slop need minus 1
*/  
forall(t in months, j in 1..t) 
  P[j][t] == 1 => stockplb[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {(- 1 + (sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 0} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock[t]/std_matrix[j][t]); 
}