/*********************************************
 * OPL 12.7.0.0 Model
 * Author: Administrator
 * Creation Date: 2018年11月13日 at 上午11:26:31
 * Description: implement the joint MILP model
 *********************************************/

 //parameters
int nbmonths=...;
range months=1..nbmonths;
float fc=...;
float v = ...;
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
dvar float stock_S[0..nbmonths];
dvar float+ stockhlb_S[0..nbmonths];
dvar float+ stockplb_S[0..nbmonths];
dvar boolean purchase_S[months];
dvar boolean P_S[months][months];

dvar float stock_s[0..nbmonths];
dvar float+ stockhlb_s[0..nbmonths];
dvar float+ stockplb_s[0..nbmonths];
dvar boolean purchase_s[months];
dvar boolean P_s[months][months];

dvar float+ C;

//float mean_matrix[i in months, j in months] = sum(m in i..j) meandemand[m];
float std_matrix[i in months, j in months] = sqrt(sum(m in i..j) pow(std_demand[m],2));

//objective function
minimize sum(t in months)( fc*purchase_S[t]+h*stockhlb_S[t]+p*stockplb_S[t]) + sum(t in 2..nbmonths)( fc*purchase_s[t]+h*stockhlb_s[t]+p*stockplb_s[t])
+ v*stock_S[nbmonths]+ v*stock_s[nbmonths];

// expressions
dexpr float Cx = sum(t in months)( fc*purchase_S[t]+h*stockhlb_S[t]+p*stockplb_S[t]) - v*stock_S[0] + v*stock_S[nbmonths];
dexpr float Gy = h*stockhlb_s[1]+p*stockplb_s[1] + sum(t in 2..nbmonths)( fc*purchase_s[t]+h*stockhlb_s[t]+p*stockplb_s[t])- v*stock_s[0] + v*stock_s[nbmonths];
dexpr float partialGy = sum(t in 2..nbmonths)( fc*purchase_s[t]+h*stockhlb_s[t]+p*stockplb_s[t]);

//constraints
subject to{
 stockhlb_S[0]==maxl(stock_S[0],0);
 stockplb_S[0]==maxl(-stock_S[0],0);
 stockhlb_s[0]==maxl(stock_s[0],0);
 stockplb_s[0]==maxl(-stock_s[0],0);
 
 purchase_S[1] == 1; // x_S[1] = 1
 purchase_s[1] == 0; // x_s[1] = 0
 
 stock_S[0] == stock_S[1] + meandemand[1];
 stock_s[0] <= stock_S[1] + meandemand[1];
 
 
 C ==  - v*stock_S[0] + sum(t in months)( fc*purchase_S[t]+h*stockhlb_S[t]+p*stockplb_S[t]) + v*stock_S[nbmonths] +v*sum(t in months)meandemand[t]; 
 C == - v*stock_s[0] + sum(t in months)( fc*purchase_s[t]+h*stockhlb_s[t]+p*stockplb_s[t]) + v*stock_s[nbmonths] +v*sum(t in months)meandemand[t];

 forall(t in months) {
   		purchase_S[t] == 0 => stock_S[t]+meandemand[t]-stock_S[t-1] == 0;  
   		purchase_s[t] == 0 => stock_s[t]+meandemand[t]-stock_s[t-1] == 0;  
 }
   		   
forall(t in months){
  		stock_S[t]+meandemand[t]-stock_S[t-1]>=0;
  		stock_s[t]+meandemand[t]-stock_s[t-1]>=0;
  }  		
  
forall (t in months){
  		sum(j in 1..t)P_S[j][t]==1;
  		sum(j in 1..t)P_s[j][t]==1;
  }  		
 
forall (t in months, j in 1..t){
  		P_S[j][t]>=purchase_S[j]-sum(k in j+1..t)purchase_S[k];
  		P_s[j][t]>=purchase_s[j]-sum(k in j+1..t)purchase_s[k];
  }  		
  
forall (t in months){
   		sum(k in 1..t) purchase_S[k] == 0 => P_S[1][t] == 1;  
   		sum(k in 1..t) purchase_s[k] == 0 => P_s[1][t] == 1;
  }
  
/**
* Original formulation as in (Rossi et al., 2015)    
*/ 
forall(t in months, p in partitions){
  stockhlb_S[t]>=sum(k in 1..p)prob[k]*stock_S[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P_S[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P_S[j][t]);
}  

forall(t in months, p in partitions)
	stockhlb_s[t]>=sum(k in 1..p)prob[k]*stock_s[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P_s[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P_s[j][t]);
  
forall(t in months) {
    stockhlb_S[t] >= (sum(j in 1..t) error*std_matrix[j][t]*P_S[j][t]);
    stockhlb_s[t] >= (sum(j in 1..t) error*std_matrix[j][t]*P_s[j][t]);
  }

    
/**
* Original formulation as in (Rossi et al., 2015)    
*/ 
forall(t in months, p in partitions){
  stockplb_S[t]>=-stock_S[t]+sum(k in 1..p)prob[k]*stock_S[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P_S[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P_S[j][t]); 
}  

forall(t in months, p in partitions)
  stockplb_s[t]>=-stock_s[t]+sum(k in 1..p)prob[k]*stock_s[t]-sum(j in 1..t)(sum(k in 1..p)prob[k]*means[k]*std_matrix[j][t]*P_s[j][t]) + (sum(j in 1..t) error*std_matrix[j][t]*P_s[j][t]);
  
 
forall(t in months) {
    stockplb_S[t] >= - stock_S[t] + (sum(j in 1..t) error*std_matrix[j][t]*P_S[j][t]);
    stockplb_s[t] >= - stock_s[t] + (sum(j in 1..t) error*std_matrix[j][t]*P_s[j][t]);
}

/**
* Piecewise formulation    
*/
forall(t in months, j in 1..t) 
  P_s[j][t] == 1 => stockhlb_s[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {((sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 1} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock_s[t]/std_matrix[j][t]);   

forall(t in months, j in 1..t) 
  P_s[j][t] == 1 => stockplb_s[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {(- 1 + (sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 0} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock_s[t]/std_matrix[j][t]); 
  
 forall(t in months, j in 1..t) 
  P_S[j][t] == 1 => stockhlb_S[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {((sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 1} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock_S[t]/std_matrix[j][t]);   

forall(t in months, j in 1..t) 
  P_S[j][t] == 1 => stockplb_S[t]/std_matrix[j][t] == 
  piecewise(i in partitions) {(- 1 + (sum(k in 1..i) prob[k]) - prob[i]) -> means[i]; 0} 
  (0, error - sum(k in 1..ftoi(round(nbpartitions/2)))(prob[k]*means[k])) 
  (stock_S[t]/std_matrix[j][t]); 
  
}

execute {
  if (cplex.getCplexStatus() == 1) {
    writeln('Cx is ', Cx);
    writeln('Gy is ', Gy);
    writeln('partial Gy is ', partialGy)
  } 
}