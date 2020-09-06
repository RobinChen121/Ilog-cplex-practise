/*********************************************
 * OPL 12.9.0.0 Model
 * Author: chen
 * Creation Date: 2020年9月6日 at 下午1:44:57
 *********************************************/

 // parameters
int items = ...;  // 每个阶段的情景数
float iniWealth = ...;
int periods = ...;
int scenarioNum = ...;
float targetWealth = ...;
float interestRate = ...; //value of investing beyond target value
float penaltyRate = ...; // penalty for not meeting target
int scenarioLinks[1..periods][1..scenarioNum][1..scenarioNum] = ...; // 0-1 array shows which scenarios are combined at period t    
float scenarioProb[1..scenarioNum] = ...;
float returns[1..periods][1..items][1..scenarioNum] = ...; // return amount in item i for scenario s in period t

// decision variables
dvar float+ x[1..periods][1..items][1..scenarioNum]; // investment amount in item i for scenario s in period t
dvar float+ y[1..scenarioNum]; // amount above final target in scenario s
dvar float+ w[1..scenarioNum]; // amount below final target in scenario s

// objective
maximize sum (i in 1..scenarioNum) scenarioProb[i]*(interestRate*y[i] - penaltyRate*w[i]);
                                                                                           
// constraints
subject to{
	iniWealthEqual:
	forall (s in 1..scenarioNum)
		sum (i in 1..items) x[1][i][s] == iniWealth;
		
	forall (s in 1..scenarioNum)
	  sum (i in 1..items) returns[periods][i][s]*x[periods][i][s] == y[s] - w[s] + targetWealth;
		
	balance: 
	forall (t in 2..periods)
	  forall (s in 1..scenarioNum)
	    sum (i in 1..items) returns[t-1][i][s]*x[t - 1][i][s] == sum (i in 1..items) x[t][i][s];
	    
	nonanticipativity:
	forall (t in 1..periods)
		forall (i in 1..items)
	      forall (s in 1..scenarioNum)  
	        (sum(s1 in 1..scenarioNum) scenarioLinks[t][s1][s]*scenarioProb[s1]*x[t][i][s1]) == 
	          (sum(s1 in 1..scenarioNum) scenarioLinks[t][s1][s]*scenarioProb[s1])*x[t][i][s]; 
}

 float y1[1..8] = [24.8, 8.87, 1.4286, 0.0, 1.4286, 0.0, 0.0, 0.0];
 float w1[1..8] = [0, 0, 0, 0.0, 0, 0.0, 0.0, 12.16];
 float x1[1..3][1..2][1..8]
  = [[[41.5, 41.5, 41.5, 41.5, 41.5, 41.5, 41.5, 41.5], [13.5, 13.5 ,13.5, 13.5, 13.5, 13.5, 13.5, 13.5]],
 	[[65.1, 65.1, 65.1, 65.1, 36.7, 36.7, 36.7, 36.7], [2.17, 2.17, 2.17, 2.17, 22.4, 22.4, 22.4, 22.4]],
	[[83.84, 83.84, 0, 0, 0, 0, 64, 64], [0, 0, 71.43, 71.43, 71.43, 71.43, 0, 0]]];
	
 
// 调试时用 
float temp = 0;
float temp1 = 0; 
float temp2 = 0;
float temp3 = 0;
range sceRange = 1..8;  // 必须用range 才能使用下面的 for 循环
range iRange = 1..2;
range tRange = 1..3;
range tRange2 = 2..3;

execute DISPLAY {  
	for (var i in sceRange) // objective
		temp += scenarioProb[i]*(interestRate*y1[i] - penaltyRate*w1[i]);		
    writeln(temp);	
    
    for (var s in sceRange) // constraint1
    {
    	temp = 0;
    	for (var i in iRange)
			temp += x1[1][i][s];	
		writeln(temp - iniWealth);		
 	}
    
    for (var s in sceRange) // constraint2
    {
    	temp = 0;
    	for (var i in iRange)
			temp += returns[periods][i][s]*x1[periods][i][s];	
		writeln(temp - y1[s] + w1[s] - targetWealth);		
 	}	
 	
 	for (var t in tRange2)	// constraint3   
    	for (var s in sceRange){ 
    		temp1 = 0; temp2 =0;	
			for (var i in iRange){
				temp = x1[t - 1][i][s]; 
				temp3 = returns[t-1][i][s];						
			 	temp1 += returns[t-1][i][s]*x1[t - 1][i][s];
			 	temp2 += x1[t][i][s];			
			}
			writeln(temp1 - temp2);
  		}				
   	
   	for (var t in tRange)	// constraint4   
    	for (var i in iRange)
			for (var s in sceRange){
				temp1 = 0; temp2 = 0;			
				for (var s1 in sceRange){
					temp1 += scenarioLinks[t][s1][s]*scenarioProb[s1]*x1[t][i][s1];
					temp2 += scenarioLinks[t][s1][s]*scenarioProb[s1];
  		 	}
  		 	temp2 = temp2*x1[t][i][s];
  		 	writeln(temp1 - temp2);		
    	}									
}